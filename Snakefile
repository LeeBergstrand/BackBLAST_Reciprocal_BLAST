# Snakefile rules for BackBLAST pipeline
# Copyright Lee Bergstrand and Jackson M. Tsuji, 2018
from snakemake.utils import logger, min_version, update_config

# Specify the minimum snakemake version allowable
min_version("5.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")


rule all:
    input:
        "generate_heatmap/BackBLAST_heatmap.pdf"


# Runs reciprocal BLAST for each subject genome against the target genes in the query genome
rule run_reciprocal_blast:
    input:
        lambda wildcards: config["subjects"][wildcards.subject]
    output:
        "reciprocal_blast/{subject}.csv"
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logs/reciprocal_blast/{subject}.log"
    benchmark:
        "benchmarks/{subject}.reciprocal_blast.benchmark.txt"
    threads: 1
    params:
        query_genes = config.get("query_genes"),
        query_genome_orfs = config.get("query_genome_orfs"),
        eval = config.get("e_value_cutoff"),
        pident = config.get("minimum_percent_identity")
    shell:
        "BackBLAST.py --gene_cluster {params.query_genes} --query_proteome {params.query_genome_orfs} --subject_proteome {input} "
            "--e_value {params.eval} --min_ident {params.pident} --output_file {output} > {log} 2>&1"


# Removes duplicate BLAST hits for each BLAST table
rule remove_duplicates:
    input:
        "reciprocal_blast/{subject}.csv"
    output:
        "remove_duplicates/{subject}.csv"
    log:
        "logs/remove_duplicates/{subject}.log"
    benchmark:
        "benchmarks/{subject}.remove_duplicates.benchmark.txt"
    threads: 1
    shell:
        "RemoveDuplicates.sh {input} > {output} 2> {log}"


# If BLAST CSV is empty, creates a blank BLAST table
rule create_blank_results:
    input:
        "remove_duplicates/{subject}.csv"
    output:
        "fix_blank_results/{subject}.csv"
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logs/create_blank_results/{subject}.log"
    benchmark:
        "benchmarks/{subject}.create_blank_results.txt"
    params:
        query_genes=config.get("query_genes")
    shell:
        "CreateBlankResults.py {input} {params.query_genes} {output} > {log} 2>&1"


# Combine the BLAST tables into a single table, and add a column for sample ID
rule combine_blast_tables:
    input:
        blast_tables=expand("fix_blank_results/{subject}.csv", subject=config.get("subjects"))
    output:
        "combine_blast_tables/blast_tables_combined.csv"
    conda:
        "envs/R_viz.yaml"
    log:
        "logs/combine_blast_tables/combine_blast_tables.log"
    benchmark:
        "benchmarks/combine_blast_tables.txt"
    shell:
        "CombineBlastTables.R {input} {output} 2> {log}"


# Generate the final heatmap
rule generate_heatmap:
    input:
        "combine_blast_tables/blast_tables_combined.csv"
    output:
        "generate_heatmap/BackBLAST_heatmap.pdf"
    conda:
        "envs/R_viz.yaml"
    log:
        "logs/generate_heatmap/generate_heatmap.log"
    benchmark:
        "benchmarks/generate_heatmap.txt"
    params:
        tree_file=config.get("phylogenetic_tree_newick"),
        tree_metadata=config.get("tree_metadata_tsv"),
        tree_decorator_colname=config.get("tree_colouring_column_name"),
        plot_custom_names=config.get("plotting_name_column"),
        gene_naming=config.get("gene_naming_tsv"),
        bootstrap_cutoff=config.get("bootstrap_cutoff"),
        root_name=config.get("root_name")
    shell:
        "generate_BackBLAST_heatmap.R -i {params.tree_file} -j {input} -o {output} -m {params.tree_metadata} "
            "-d {params.tree_decorator_colname} -n {params.plot_custom_names} -g {params.gene_naming} "
            "-b {params.bootstrap_cutoff} -r {params.root_name} 2> {log}"
