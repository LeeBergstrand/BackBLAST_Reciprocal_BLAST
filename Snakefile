configfile: "config.yaml" # TODO - specify this during running the script instead of hard-coding

# Specify the minimum snakemake version allowable
min_version("5.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")

rule all:
    input:
        "report.html"

# Runs reciprocal BLAST for each subject genome against the target genes in the query genome
rule run_reciprocal_blast
    input:
        lambda wildcards: config["subjects"][wildcards.subject]
    output:
        "reciprocal_blast/{subject}.tsv"
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logs/reciprocal_blast/{subject}.log"
    benchmark:
        "benchmarks/{subject}.reciprocal_blast.benchmark.txt"
    threads: 1 # TODO - can it support more than one thread? Also, should I add a memory setting?
    params:
        query_genes=config.get("query_genes")
        query_genome_orfs=config.get("query_genome_orfs")
        eval = config.get("e_value_cutoff", 0.000001)
        pident = config.get("minimum_percent_identity", 25)
    shell:
        # TODO - make sure the flags match the real flags.
        "BackBLAST.py --query_genes {params.query_genes} --query_orfs {params.query_genome_orfs} --subject {input} "
            "--eval {params.eval} --pident {params.pident} --threads {threads} > {output} 2> {log}"

# Removes duplicate BLAST hits for each subject??
rule remove_duplicates
    input:
        "reciprocal_blast/{subject}.tsv"
    output:
        "remove_duplicates/{subject}.tsv"
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logs/remove_duplicates/{subject}.log"
    benchmark:
        "benchmarks/{subject}.remove_duplicates.benchmark.txt"
    threads:
        config.get("threads", 1)
    params:
        query_genes = config.get("query_genes")
        query_genome_orfs = config.get("query_genome_orfs")
    shell:
    # TODO - make sure the flags match the real structure.
    "Visualization/RemoveDuplicates.sh {input} > {output} 2> {log}"
