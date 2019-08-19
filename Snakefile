# Snakefile rules for BackBLAST pipeline
# Copyright Lee Bergstrand and Jackson M. Tsuji, 2018
import os
from snakemake.utils import logger, min_version, update_config

# Specify the minimum snakemake version allowable
min_version("5.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")

rule all:
    input:
        "logging/checkpoints/finished_blast",
        "logging/checkpoints/finished_heatmap"

# Run reciprocal BLAST for each subject genome against the target genes in the query genome
rule run_reciprocal_blast:
    input:
        lambda wildcards: config["subjects"][wildcards.subject]
    output:
        temp("blast/intermediate/reciprocal_blast/{subject}.csv")
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logging/logs/blast/intermediate/reciprocal_blast/{subject}.log"
    benchmark:
        "logging/benchmarks/{subject}.reciprocal_blast.benchmark.txt"
    threads: 1
    params:
        query_genes = config.get("query_genes"),
        query_genome_orfs = config.get("query_genome_orfs"),
        eval = config.get("e_value_cutoff"),
        pident = config.get("minimum_percent_identity")
    shell:
        "BackBLAST.py --gene_cluster {params.query_genes} --query_proteome {params.query_genome_orfs} --subject_proteome {input} "
            "--e_value {params.eval} --min_ident {params.pident} --output_file {output} > {log} 2>&1"

# Remove duplicate BLAST hits for each BLAST table
rule remove_duplicates:
    input:
        "blast/intermediate/reciprocal_blast/{subject}.csv"
    output:
        temp("blast/intermediate/remove_duplicates/{subject}.csv")
    log:
        "logging/logs/blast/intermediate/remove_duplicates/{subject}.log"
    benchmark:
        "logging/benchmarks/{subject}.remove_duplicates.benchmark.txt"
    threads: 1
    shell:
        "RemoveDuplicates.sh {input} > {output} 2> {log}"

# If BLAST CSV is empty, create a blank BLAST table
rule create_blank_results:
    input:
        "blast/intermediate/remove_duplicates/{subject}.csv"
    output:
        temp("blast/intermediate/fix_blank_results/{subject}.csv")
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logging/logs/blast/intermediate/create_blank_results/{subject}.log"
    benchmark:
        "logging/benchmarks/{subject}.create_blank_results.benchmark.txt"
    params:
        query_genes=config.get("query_genes")
    shell:
        "CreateBlankResults.py -i {input} -q {params.query_genes} -o {output} > {log} 2>&1"

# Combine the BLAST tables into a single table, and add a column for sample ID
rule combine_blast_tables:
    input:
        blast_tables=expand("blast/intermediate/fix_blank_results/{subject}.csv", subject=config.get("subjects"))
    output:
        "blast/combine_blast_tables/blast_tables_combined.csv"
    conda:
        "envs/R_viz.yaml"
    log:
        "logging/logs/blast/combine_blast_tables/combine_blast_tables.log"
    benchmark:
        "logging/benchmarks/combine_blast_tables.benchmark.txt"
    shell:
        "CombineBlastTables.R {input} {output} > {log} 2>&1"

# Create a symlink of the combined BLAST table to make it easier for the user to find
rule symlink_combined_blast_table:
    input:
        "blast/combine_blast_tables/blast_tables_combined.csv"
    output:
        "blast/blast_tables_combined.csv"
    run:
        source_relpath = os.path.relpath(str(input), os.path.dirname(str(output)))
        os.symlink(source_relpath, str(output))

# Checkpoint that BLAST step is finished
rule finished_blast:
    input: "blast/blast_tables_combined.csv"
    output: touch("logging/checkpoints/finished_blast")

# Generate phylogenetic tree if desired by the user
if config.get("phylogenetic_tree_newick") == "subjects":
    # Makes a list of filepaths to the FAA genome files for use by GToTree
    # Part 1 - generate for each sample
    rule generate_phylogenetic_tree_input_individual:
        input:
            lambda wildcards: config["subjects"][wildcards.subject]
        output:
            temp("phylogeny/input/{subject}.list")
        threads: 1
        shell:
            "echo {input} > {output}"

    # Part 2 of the above rule
    # TODO - check {output} does not already exist before run start?
    # TODO - how to merge parts 1 and 2 together?
    rule generate_phylogenetic_tree_input_grouped:
        input:
            blast_tables=expand("phylogeny/input/{subject}.list", subject=config.get("subjects"))
        output:
            "phylogeny/input/input_genomes_faa.list"
        threads: 1
        shell:
            "cat {input} >> {output}"

    # Run GToTree
    # TODO - expose more params to the user
    rule generate_phylogenetic_tree:
        input:
            "phylogeny/input/input_genomes_faa.list"
        output:
            "phylogeny/iqtree_out.treefile"
        conda:
            "envs/gtotree.yaml"
        log:
            "logging/logs/phylogeny/gtotree.log"
        benchmark:
            "logging/benchmarks/gtotree.benchmark.txt"
        threads: config.get("threads", 1)
        params:
            phylogenetic_model = config.get("gtotree_phylogenetic_model", "Universal_Hug_et_al.hmm"),
            sequence_length_threshold = config.get("gtotree_sequence_length_threshold", "0.2"),
            minimum_hit_fraction = config.get("gtotree_minimum_hit_fraction", "0.5")
        shell:
            "GToTree -A {input} -H {params.phylogenetic_model} -o phylogeny/gtotree -T IQ-TREE "
                "-c {params.sequence_length_threshold} -G {params.minimum_hit_fraction} "
                "-n {threads} -j {threads} > {log} 2>&1 && "
            "ln phylogeny/gtotree/iqtree_out/iqtree_out.treefile phylogeny/iqtree_out.treefile"

# Create a fake temp file to allow plotter to run if 'NA' is selected
# TODO - this is a hack; is there a proper way to make this work?
if config.get("phylogenetic_tree_newick") == "NA":
    rule generate_fake_phylogenetic_tree:
        output: temp(touch("NA"))

# Generate the final heatmap
rule generate_heatmap:
    input:
        blast_table = "blast/combine_blast_tables/blast_tables_combined.csv",
        tree_file = "phylogeny/iqtree_out.treefile" if config.get("phylogenetic_tree_newick") == "subjects" else config.get("phylogenetic_tree_newick", "NA")
    output:
        "heatmap/BackBLAST_heatmap.pdf"
    conda:
        "envs/R_viz.yaml"
    log:
        "logging/logs/heatmap/generate_heatmap.log"
    benchmark:
        "logging/benchmarks/generate_heatmap.benchmark.txt"
    params:
        genome_metadata = config.get("genome_metadata_tsv", "NA"),
        gene_metadata = config.get("gene_metadata_tsv", "NA"),
        bootstrap_cutoff = config.get("bootstrap_cutoff", "NA"),
        root_name = config.get("root_name", "NA"),
        plot_width = config.get("plot_width_mm", 400),
        plot_height = config.get("plot_height_mm", 200)
    shell:
        "generate_BackBLAST_heatmap.R -m {params.genome_metadata} -g {params.gene_metadata} "
            "-b {params.bootstrap_cutoff} -r {params.root_name} -w {params.plot_width} -z {params.plot_height} "
            "{input.tree_file} {input.blast_table} {output} 2>&1 | tee {log} && "
        "if [[ -f Rplots.pdf ]]; then rm Rplots.pdf; fi"

# Checkpoint that the plot is done
# This is probably over-engineered, but it is helpful to mirror what is being done for the BLAST step and could help to one day separate this Snakefile into modules
rule finished_heatmap:
    input: "heatmap/BackBLAST_heatmap.pdf"
    output: touch("logging/checkpoints/finished_heatmap")

