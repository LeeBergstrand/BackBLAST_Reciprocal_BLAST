#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# Copyright: Lee H. Bergstrand (2024)
# Description: A Biopython program that takes a list of query proteins and uses local BLASTP to search
#              for highly similar proteins within a local blast database (usually a local db of a target
#              proteome). The program then BLASTPs backwards from the found subject proteins to the query
#              proteome to confirm gene orthology. Part of the BackBLAST pipeline.
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
import networkx as nx

DEFAULT_E_VALUE_CUTOFF = 1e-25
DEFAULT_MINIMUM_IDENTITY_CUTOFF = 25
DEFAULT_MINIMUM_QUERY_COVERAGE = 50


def get_blast_hight_scoring_pairs(query_gene_cluster_path, subject_proteome_file, e_value_cutoff, minimum_identity,
                                  minimum_query_coverage):
    """
    Runs BLASTP on the given query and subject FASTA files and returns collects high scoring pairs.


    :param query_gene_cluster_path: The path to the query FASTA file.
    :param subject_proteome_file: The path to the subject FASTA file.
    :param e_value_cutoff: The e-value cutoff for BLASTP.
    :param minimum_identity: The minimum sequence identity that each HSP should have.
    :param minimum_query_coverage: The minimum query coverage (qcovhsp) that each HSP should have.
    :return: The blast output as a list of lists representing HSPs and their parameters.
    """
    return filter_blast_csv(run_blastp(query_gene_cluster_path,
                                       subject_proteome_file,
                                       e_value_cutoff=e_value_cutoff),
                            minimum_identity=minimum_identity,
                            minimum_query_coverage=minimum_query_coverage)


def run_blastp(query_file_path, subject_file_path, e_value_cutoff):
    """
    Runs BLASTP on the given query and subject FASTA files.

    :param query_file_path: The path to the query FASTA file.
    :param subject_file_path: The path to the subject FASTA file.
    :param e_value_cutoff: The e-value cutoff for BLASTP.
    :return:    A csv formatted BLASTP output (query_sequence_id, subject_sequence_id, percent_identity, e-value,
                query coverage, bitscore)
    """
    blast_out = subprocess.check_output(
        ["blastp", "-query", query_file_path, "-subject", subject_file_path, "-evalue", str(e_value_cutoff),
         "-soft_masking", "true", "-seg", "yes", "-outfmt", "10 qseqid sseqid pident evalue qcovhsp bitscore"])

    # Decodes BLASTP output to UTF-8 (In Py3 check_output returns raw bytes)
    blast_out = blast_out.decode().replace(' ', '')
    return blast_out


def filter_blast_csv(raw_blast_output, minimum_identity, minimum_query_coverage):
    """
    This filters the blast results by minimum_identity and creates a 2D array of HSPs.

    :param raw_blast_output: String containing the entire output from BLAST.
    :param minimum_identity: The minimum sequence identity that each HSP should have.
    :param minimum_query_coverage: The minimum query coverage (qcovhsp) that each HSP should have.
    :return: The filter blast output as a list of lists representing HSPs and their parameters.
    """
    blast_results = csv.reader(raw_blast_output.splitlines(True))  # Reads BLAST csv rows as a csv

    filtered_blast_output = []  # Note should simply delete unwanted HSPs from current list rather than making new list
    # Rather than making a new one
    for high_scoring_pair in blast_results:
        if float(high_scoring_pair[2]) >= minimum_identity:  # Filter by minimum identity
            if float(high_scoring_pair[4]) >= minimum_query_coverage:  # Filter by minimum query coverage
                # Converts each high_scoring_pair parameter that should be a number to a number
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
        handle = open(file_path, "r")
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
        print(hit)
        subject_protein = blast_graph[hit[1]]

        # Find the top score of the best reciprocal BLAST hit
        top_back_hit_score = 0
        for back_hit_id in subject_protein:
            back_hit_score = get_edge_attribute(subject_protein, back_hit_id,
                                                'bitscore-rev', none_value=0)
            if back_hit_score >= top_back_hit_score:
                top_back_hit_score = back_hit_score

        # Check if the query is the best reciprocal BLAST hit for the subject
        delete_hit = False
        if hit[0] in subject_protein:
            # The edge weight between the subject and the query is the reciprocal BLAST score
            back_hit_to_query_score = get_edge_attribute(subject_protein, hit[0],
                                                         'bitscore-rev', none_value=0)
            if back_hit_to_query_score < top_back_hit_score:
                # If the query is not the best reciprocal BLAST hit simply delete
                # it from the filterable_forward_blast_results
                delete_hit = True
        else:
            # If the query is not a reciprocal BLAST hit simply delete it from the filterable_forward_blast_results
            delete_hit = True

        if delete_hit:
            # Delete the forward BLAST hit from filterable_forward_blast_results
            del filterable_forward_blast_results[filterable_forward_blast_results.index(hit)]

    # Return reciprocal BLAST output
    return filterable_forward_blast_results


def create_blast_graph(forward_blast_high_scoring_pairs, reverse_blast_high_scoring_pairs):
    """
    Builds a graph network of forward and reverse hits.

    :param forward_blast_high_scoring_pairs: The forward blast high scoring pairs.
    :param reverse_blast_high_scoring_pairs: The reverse blast high scoring pairs.
    :return: A graph representation of forward and reverse HSPs.
    """
    blast_graph = nx.Graph()  # Creates graph to map BLAST hits
    print(">> Creating Graph...")
    # hit[0] = query ID; hit[1] = subject ID; hit[5] = bitscore

    for hit in forward_blast_high_scoring_pairs:
        blast_graph.add_edge(hit[0], hit[1])
        blast_graph[hit[0]][hit[1]]['bitscore-fwd'] = hit[5]

    for hit in reverse_blast_high_scoring_pairs:
        blast_graph.add_edge(hit[1], hit[0])
        blast_graph[hit[1]][hit[0]]['bitscore-rev'] = hit[5]

    return blast_graph


def get_edge_attribute(vertex, neighbour_id, attribute_name, none_value=None):
    """
    Get attribute of the edge of a nx.Graph()

    :param vertex: AtlasView of a graph centered on the vertex of interest (e.g., graph[vertex]).
    :param neighbour_id: The ID of the neighbouring vertex to form the edge.
    :param attribute_name: Name of the edge attribute.
    :param none_value: What to return of the attribute does not exist for the edge.
    :return: The value of the attribute for the desired edge.
    """
    edge_attributes = vertex[neighbour_id]

    if attribute_name in edge_attributes:
        attribute_value = edge_attributes[attribute_name]
    else:
        attribute_value = none_value

    return attribute_value


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
    input_min_query_cov_cutoff = args.min_query_cov
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
    # Forward BLAST from query proteins to subject proteome and filter BLAST results by percent identity and query cov
    forward_blast_high_scoring_pairs = get_blast_hight_scoring_pairs(query_gene_cluster_path=query_gene_cluster_path,
                                                                     subject_proteome_file=subject_proteome_file,
                                                                     e_value_cutoff=input_e_value_cutoff,
                                                                     minimum_identity=input_min_ident_cutoff,
                                                                     minimum_query_coverage=input_min_query_cov_cutoff)

    if len(forward_blast_high_scoring_pairs) == 0:
        print(">> No Forward hits in subject proteome were found.")

        try:
            open(out_file, "w").close()  # Writes empty file for easier data processing.
        except IOError:
            print(">> Failed to create " + out_file)
            sys.exit(1)
        print(">> Exiting.\n\n")
        sys.exit(0)  # Aborts program. (exit(0) indicates that no error occurred)

    # Creates python dictionary with every protein FASTA sequence in the subject proteome
    subject_proteome_fasta_cache = create_fasta_cache(subject_proteome_file)

    print(">> Creating Back-Blasting Query from found subject proteins...")
    # For each top hit...
    back_blast_query_fastas = []
    for hit in forward_blast_high_scoring_pairs:
        subject_protein = hit[1]
        subject_protein_fasta = subject_proteome_fasta_cache.get(subject_protein)
        back_blast_query_fastas.append(subject_protein_fasta)  # Adds current subject to overall protein list

    complete_back_blast_query = "".join(back_blast_query_fastas)

    # Attempt to write a temporary FASTA file for the reverse BLAST to use
    temp_filename = "temp_query_" + uuid.uuid4().hex + ".faa"
    print(">> Writing backBLASTing query to temporary file " + temp_filename)
    try:
        write_file = open(temp_filename, "w")
        write_file.write(complete_back_blast_query)
        write_file.close()
    except IOError:
        print("Failed to create " + temp_filename)
        sys.exit(1)

    print(">> BLASTing backwards from subject genome to query genome.")
    # Run backwards BLAST towards query proteome and filters BLAST results by percent identity
    reverse_blast_high_scoring_pairs = get_blast_hight_scoring_pairs(query_gene_cluster_path=temp_filename,
                                                                     subject_proteome_file=query_proteome_path,
                                                                     e_value_cutoff=input_e_value_cutoff,
                                                                     minimum_identity=input_min_ident_cutoff,
                                                                     minimum_query_coverage=input_min_query_cov_cutoff)

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

    parser.add_argument('-e', '--e_value', metavar='E-VALUE', default=DEFAULT_E_VALUE_CUTOFF, type=str,
                        help='''The Expect value (E) cutoff for removing high scoring pairs. ''' +
                             '''The smaller this number is, the stricter the BLAST search.''')

    parser.add_argument('-i', '--min_ident', metavar='IDENT', default=DEFAULT_MINIMUM_IDENTITY_CUTOFF, type=float,
                        help='''The minimum sequence identify cutoff for removing high scoring pairs.''' +
                             '''The larger this number is, the stricter the BLAST search.''')

    parser.add_argument('-c', '--min_query_cov', metavar='QCOV', default=DEFAULT_MINIMUM_QUERY_COVERAGE, type=float,
                        help='''The minimum percent query coverage cutoff for removing high scoring pairs.''' +
                             '''The larger this number is, the stricter the BLAST search.''')

    parser.add_argument('-o', '--output_file', metavar='OUTPUT', required=True,
                        help='''The path to write CSV-format BLAST results to.''')

    cli_args = parser.parse_args()
    main(cli_args)
