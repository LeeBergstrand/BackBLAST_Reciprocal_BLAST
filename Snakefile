configfile: "config.yaml"

# Specify the minimum snakemake version allowable
min_version("5.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")

rule all:
    input:
        "report.html"

rule run_reciprocal_blast
    input:
        genomes = config.get("genome_ORFs_tsv")
        genes = config.get("gene_targets_tsv")
    output:
        "reciprocal_blast/reciprocal_blast_results.tsv"
    conda:
        "envs/reciprocal_blast.yaml"
    log:
        "logs/reciprocal_blast.log"
    benchmark:
        "benchmarks/reciprocal_blast.benchmark.txt"
    threads:
        config.get("threads", 1)
    params:
        eval = config.get("e_value_cutoff", 0.000001)
        pident = config.get("minimum_percent_identity", 25)
    shell:
        # TODO - fix this. I just threw in a template command for now, but the required inputs don't even match those that I defined above! I have the command as I'd see it ideally as the second line below, with the real command above.
        "BackBLAST.py <queryGeneList.faa> <queryBLASTDB.faa> <subjectBLASTDB.faa> > {output} 2> {log}"
        "BackBLAST.py --genes {input.genes} --orfs {input.genomes} --eval {params.eval} --pident {params.pident} --threads {threads} > {output} 2> {log}"

