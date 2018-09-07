configfile: "config.yaml" # TODO - specify this during running the script instead of hard-coding

# Specify the minimum snakemake version allowable
min_version("5.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")

rule all:
    input:
        "report.html"

rule run_first_blast
    input:
        lambda wildcards: config["subjects"][wildcards.subject]
    output:
        "blast_1/{subject}.tsv"
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logs/blast_1/{subject}.log"
    benchmark:
        "benchmarks/{subject}.blast_1.benchmark.txt"
    threads:
        config.get("threads", 1)
    params:
        query_genes=config.get("query_genes")
        query_genome_orfs=config.get("query_genome_orfs")
        eval = config.get("e_value_cutoff", 0.000001)
        pident = config.get("minimum_percent_identity", 25)
    shell:
        # TODO - make sure the flags match the real flags.
        "BackBLAST.py --query_genes {params.query_genes} --query_orfs {params.query_genome_orfs} --subject {input} "
            "--eval {params.eval} --pident {params.pident} --threads {threads} > {output} 2> {log}"

