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
        "reciprocal_blast/{subject}.csv"
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
        "BackBLAST.py --gene_cluster {params.query_genes} --query_proteome {params.query_genome_orfs} --subject_proteome {input} "
            "--eval {params.eval} --pident {params.pident} --threads {threads} > {output} 2> {log}"

# Removes duplicate BLAST hits for each BLAST table
rule remove_duplicates
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
    # TODO - make sure the flags match the real structure.
    "Visualization/RemoveDuplicates.sh {input} > {output} 2> {log}"

# TODO - this is only run on SOME outputs, correct? How to construct an IF statement here?
rule create_blank_results
    input:
        "remove_duplicates/{subject}.csv"
    output:
        "fix_blank_results/{subject}.csv"
    conda:
        "envs/reciprocal_blast.yaml"
    params:
        query_genes = config.get("query_genes")
    shell:
    # TODO - original script needs 'output'!
        "CreateBlankResults.py {input} {output} {params.query_genes}"