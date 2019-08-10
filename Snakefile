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
        "CreateBlankResults.py -i {input} -q {params.query_genes} -o {output} > {log} 2>&1"


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
        "CombineBlastTables.R {input} {output} > {log} 2>&1"


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
            "logs/phylogeny/gtotree.log"
        benchmark:
            "benchmarks/gtotree.txt"
        threads: config.get("threads", 1)
        params:
            phylogenetic_model = "Universal_Hug_et_al.hmm"
        shell:
            "GToTree -A {input} -H {params.phylogenetic_model} -o phylogeny/gtotree -T IQ-TREE -c 0.2 -G 0.5 "
                "-n {threads} -j {threads} > {log} 2>&1 && "
            "ln phylogeny/gtotree/iqtree_out/iqtree_out.treefile phylogeny/iqtree_out.treefile"


# Generate the final heatmap
rule generate_heatmap:
    input:
        blast_table = "combine_blast_tables/blast_tables_combined.csv",
        tree_file = "phylogeny/iqtree_out.treefile" if config.get("phylogenetic_tree_newick") == "subjects" else config.get("phylogenetic_tree_newick")
    output:
        "generate_heatmap/BackBLAST_heatmap.pdf"
    conda:
        "envs/R_viz.yaml"
    log:
        "logs/generate_heatmap/generate_heatmap.log"
    benchmark:
        "benchmarks/generate_heatmap.txt"
    params:
        genome_metadata = config.get("genome_metadata_tsv", "NA"),
        gene_metadata = config.get("gene_metadata_tsv", "NA"),
        bootstrap_cutoff = config.get("bootstrap_cutoff", "NA"),
        root_name = config.get("root_name", "NA")
    shell:
        "generate_BackBLAST_heatmap.R -m {params.genome_metadata} -g {params.gene_metadata} "
            "-b {params.bootstrap_cutoff} -r {params.root_name} {input.tree_file} {input.blast_table} {output} 2>&1 | tee {log}"

