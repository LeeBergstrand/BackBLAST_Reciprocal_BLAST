#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# Copyright: Lee H. Bergstrand (2019)
# Description: A Biopython program that takes a list of query proteins and uses local BLASTp to search
#              for highly similar proteins within a local blast database (usually a local db of a target
#              proteome). The program then BLASTps backwards from the found subject proteins to the query
#              proteome to confirm gene orthology.
#             
# Requirements: - This script requires BLAST+ 2.2.9 or later.
#               - All operations are done with protein sequences.
#               - All query proteins should be from sequenced genomes in order to facilitate backwards BLAST.
#
# ===========================================================================================================

import argparse
# Imports & Setup:
import csv
import subprocess
import sys
import uuid
import os

from Bio import SeqIO

from Graph import Graph

DEFAULT_E_VALUE_CUTOFF = 1e-25
DEFAULT_MINIMUM_IDENTITY_CUTOFF = 25


def get_blast_hight_scoring_pairs(query_gene_cluster_path, subject_proteome_file, e_value_cutoff, minimum_identity):
    """
    Runs BLASTp on the given query and subject FASTA files and returns collects high scoring pairs.


    :param query_gene_cluster_path: The path to the query FASTA file.
    :param subject_proteome_file: The path to the subject FASTA file.
    :param e_value_cutoff: The e-value cutoff for BLASTp.
    :param minimum_identity: The minimum sequence identity that each HSP should have.
    :return: The blast output as a list of lists representing HSPs and their parameters.
    """
    return filter_blast_csv(run_blastp(query_gene_cluster_path,
                                       subject_proteome_file,
                                       e_value_cutoff=e_value_cutoff),
                            minimum_identity=minimum_identity)


def run_blastp(query_file_path, subject_file_path, e_value_cutoff):
    """
    Runs BLASTp on the given query and subject FASTA files.

    :param query_file_path: The path to the query FASTA file..
    :param subject_file_path: The path to the subject FASTA file.
    :param e_value_cutoff: The e-value cutoff for BLASTp.
    :return:    A csv formatted BLASTp output (query_sequence_id, subject_sequence_id, percent_identity, e-value,
                query coverage, bitscore)
    """
    blast_out = subprocess.check_output(
        ["blastp", "-query", query_file_path, "-subject", subject_file_path, "-evalue", str(e_value_cutoff),
         "-soft_masking", "true", "-seg", "yes", "-outfmt", "10 qseqid sseqid pident evalue qcovhsp bitscore"])

    # Decodes BLASTp output to UTF-8 (In Py3 check_output returns raw bytes)
    blast_out = blast_out.decode().replace(' ', '')
    return blast_out


def filter_blast_csv(raw_blast_output, minimum_identity):
    """
    This filters the blast results by minimum_identity and creates a 2D array of HSPs.

    :param raw_blast_output: String containing the entire output from BLAST.
    :param minimum_identity: The minimum sequence identity that each HSP should have.
    :return: The filter blast output as a list of lists representing HSPs and their parameters.
    """
    blast_results = csv.reader(raw_blast_output.splitlines(True))  # Reads BLAST csv rows as a csv.

    filtered_blast_output = []  # Note should simply delete unwanted HSPs from current list rather than making new list.
    # Rather than making a new one.
    for high_scoring_pair in blast_results:
        if float(high_scoring_pair[2]) >= minimum_identity:  # Filter by minimum identity.
            # Converts each high_scoring_pair parameter that should be a number to a number.
            high_scoring_pair[2] = float(high_scoring_pair[2])
            high_scoring_pair[3] = float(high_scoring_pair[3])
            high_scoring_pair[4] = float(high_scoring_pair[4])
            high_scoring_pair[5] = float(high_scoring_pair[5])
            filtered_blast_output.append(high_scoring_pair)

    return filtered_blast_output


def create_fasta_cache(file_path):
    """
    Creates a cache of FASTA sequences from a FASTA file.

    :param file_path: A path to a FASTA file.
    :return: A dictionary of FASTA sequences keyed by each sequence's ID.
    """
    proteome_hash = dict()
    try:
        handle = open(file_path, "rU")
        for record in SeqIO.parse(handle, "fasta"):
            proteome_hash.update({record.id: record.format("fasta")})
        handle.close()
    except IOError:
        print("Failed to open " + file_path)
        sys.exit(1)

    return proteome_hash


def filter_forward_pairs_by_reverse_pairs(forward_blast_high_scoring_pairs, reverse_blast_high_scoring_pairs):
    """
    Takes the forward blast results and uses the reverse blast results to filter out non-orthologus HSPs.

    Uses the algorithm found here:

    https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST/wiki/BackBLAST-Algorithm

    :param forward_blast_high_scoring_pairs: The forward blast high scoring pairs.
    :param reverse_blast_high_scoring_pairs: The reverse blast high scoring pairs.
    :return: The forward blast results with non-ortholgous HSPs removed.
    """
    blast_graph = create_blast_graph(forward_blast_high_scoring_pairs, reverse_blast_high_scoring_pairs)

    filterable_forward_blast_results = list(forward_blast_high_scoring_pairs)
    print(">> Checking if forward hit subjects have better reciprocal hits than query.")
    for hit in forward_blast_high_scoring_pairs:
        query_protein = blast_graph.getVertex(hit[0])
        subject_protein = blast_graph.getVertex(hit[1])

        top_back_hit_score = 0
        # Find the top score of the best reciprocal BLAST hit.
        for backHit in subject_protein.getConnections():
            back_hit_score = subject_protein.getWeight(
                backHit)  # The edge weight between the subject and its reciprocal BLAST hit is the BLAST score.
            if back_hit_score >= top_back_hit_score:
                top_back_hit_score = back_hit_score

        # Check if the query is the best reciprocal BLAST hit for the subject.
        delete_hit = False
        if query_protein in subject_protein.getConnections():
            # The edge weight between the subject and the query is the reciprocal BLAST score.
            back_hit_to_query_score = subject_protein.getWeight(query_protein)

            if back_hit_to_query_score < top_back_hit_score:
                # If the query is not the best reciprocal BLAST hit simply delete
                # it from the filterable_forward_blast_results.
                delete_hit = True
        else:
            # If the query is not a reciprocal BLAST hit simply delete it from the filterable_forward_blast_results.
            delete_hit = True

        if delete_hit:
            # Delete the forward BLAST hit from filterable_forward_blast_results.
            del filterable_forward_blast_results[filterable_forward_blast_results.index(hit)]

            # Attempts to write reciprocal BLAST output to file.
    return filterable_forward_blast_results


def create_blast_graph(forward_blast_high_scoring_pairs, reverse_blast_high_scoring_pairs):
    """
    Builds a graph network of forward and reverse hits.

    :param forward_blast_high_scoring_pairs: The forward blast high scoring pairs.
    :param reverse_blast_high_scoring_pairs: The reverse blast high scoring pairs.
    :return: A graph representation of forward and reverse HSPs.
    """
    blast_graph = Graph()  # Creates graph to map BLAST hits.
    print(">> Creating Graph...")

    for hit in forward_blast_high_scoring_pairs:
        blast_graph.addEdge(hit[0], hit[1], hit[5])
    for hit in reverse_blast_high_scoring_pairs:
        blast_graph.addEdge(hit[0], hit[1], hit[5])

    return blast_graph


def main(args):
    """
    The starting point of BackBLAST.

    :param args: The CLI arguments.
    """
    query_gene_cluster_path = args.gene_cluster
    query_proteome_path = args.query_proteome
    subject_proteome_file = args.subject_proteome
    input_e_value_cutoff = args.e_value
    input_min_ident_cutoff = args.min_ident
    out_file = args.output_file

    print("Opening " + subject_proteome_file + "...")

    # File extension checks
    if not query_gene_cluster_path.endswith(".faa"):
        print("[Warning] " + query_gene_cluster_path + " may not be a amino acid FASTA file!")
    if not query_proteome_path.endswith(".faa"):
        print("[Warning] " + query_proteome_path + " may not be a amino acid FASTA file!")
    if not subject_proteome_file.endswith(".faa"):
        print("[Warning] " + subject_proteome_file + " may not be a amino acid FASTA file!")


    print(">> Forward Blasting to subject proteome...")
    # Forward BLASTs from query proteins to subject proteome and filters BLAST results by percent identity.
    forward_blast_high_scoring_pairs = get_blast_hight_scoring_pairs(query_gene_cluster_path=query_gene_cluster_path,
                                                                      subject_proteome_file=subject_proteome_file,
                                                                      e_value_cutoff=input_e_value_cutoff,
                                                                      minimum_identity=input_min_ident_cutoff)

    if len(forward_blast_high_scoring_pairs) == 0:
        print(">> No Forward hits in subject proteome were found.")

        try:
            open(out_file, "w").close()  # Writes empty file for easier data processing.
        except IOError:
            print(">> Failed to create " + out_file)
            sys.exit(1)
        print(">> Exiting.\n\n")
        sys.exit(0)  # Aborts program. (exit(0) indicates that no error occurred)

    # Creates python dictionary with every protein FASTA sequence in the subject proteome.
    subject_proteome_fasta_cache = create_fasta_cache(subject_proteome_file)

    print(">> Creating Back-Blasting Query from found subject proteins...")
    # For each top Hit...
    back_blast_query_fastas = []
    for hit in forward_blast_high_scoring_pairs:
        subject_protein = hit[1]
        subject_protein_fasta = subject_proteome_fasta_cache.get(subject_protein)
        back_blast_query_fastas.append(subject_protein_fasta)  # Adds current subject to overall protein list.

    complete_back_blast_query = "".join(back_blast_query_fastas)

    # Attempt to write a temporary FASTA file for the reverse BLAST to use.
    try:
        temp_filename = "temp_query_" + uuid.uuid4().hex + ".faa"
        print(">> Writing Back-Blasting Query to temporary file " + temp_filename)
        write_file = open(temp_filename, "w")
        write_file.write(complete_back_blast_query)
        write_file.close()
    except IOError:
        print("Failed to create " + temp_filename)
        sys.exit(1)

    print(">> Blasting backwards from subject genome to query genome.")
    # Run backwards BLAST towards query proteome and filters BLAST results by percent identity.
    reverse_blast_high_scoring_pairs = get_blast_hight_scoring_pairs(query_gene_cluster_path=temp_filename,
                                                                      subject_proteome_file=query_proteome_path,
                                                                      e_value_cutoff=input_e_value_cutoff,
                                                                      minimum_identity=input_min_ident_cutoff)

    filterable_forward_blast_results = filter_forward_pairs_by_reverse_pairs(forward_blast_high_scoring_pairs,
                                                                             reverse_blast_high_scoring_pairs)
    try:
        write_file = open(out_file, "w")
        writer = csv.writer(write_file)
        print(">> Output file created.")
        print(">> Writing Data...")
        for row in filterable_forward_blast_results:
            writer.writerow(row)
        write_file.close()
        os.remove(temp_filename)
    except IOError:
        print(">> Failed to create " + out_file)
        sys.exit(1)
    print(">> Done\n")


if __name__ == '__main__':
    """Command Line Interface Options"""

    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--gene_cluster', metavar='FASTA', required=True,
                        help='''The path to the protein FASTA file of the gene cluster to be used as a query.''')

    parser.add_argument('-r', '--query_proteome', metavar='FASTA', required=True,
                        help='''The path to a FASTA file containing all proteins from query organism.''')

    parser.add_argument('-s', '--subject_proteome', metavar='FASTA', required=True,
                        help='''The path to a FASTA file containing all proteins from subject organism.''')

    parser.add_argument('-e', '--e_value', metavar='E-VALUE', default=DEFAULT_E_VALUE_CUTOFF, type=float,
                        help='''The Expect value (E) cutoff for removing high scoring pairs. ''' +
                             '''The smaller this number is, the stricter the BLAST search.''')

    parser.add_argument('-i', '--min_ident', metavar='IDENT', default=DEFAULT_MINIMUM_IDENTITY_CUTOFF, type=float,
                        help='''The minimum sequence identify cutoff for removing high scoring pairs.''' +
                             '''The larger this number is, the stricter the BLAST search.''')

    parser.add_argument('-o', '--output_file', metavar='OUTPUT', required=True,
                        help='''The path to write CSV-format BLAST results to.''')

    cli_args = parser.parse_args()
    main(cli_args)
