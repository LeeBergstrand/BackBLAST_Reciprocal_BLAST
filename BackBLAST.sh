#!/usr/bin/env bash
set -euo pipefail
# BackBLAST.sh
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019
# This is the entry script for running BackBLAST

# TODO - check for this and add a proper help message and/or parser
run_mode=$1 # setup, run, or auto

# Function to add the user-provided samples to the end of the template config file
function add_samples_to_config_file {
    # User-provided inputs
    output_config_filepath=$1
    query_faa_directory=$2
    query_reference_faa_directory=$3
    subject_faa_directory=$4

    # Query faa files
    printf "query_genes:" >> ${output_config_filepath}
    local query_faa_files=($(find ${query_faa_directory} -maxdepth 1 -type f -name "*.faa" | sort -h | xargs realpath))
    # TODO - confirm there is only one
    printf "${query_faa_files}\n" >> ${output_config_filepath}
    #for query_faa_file in ${query_faa_files[@]}; do

    #    query_faa_basename=${query_faa_file%gene_targets.faa}
    #    query_faa_basename=${query_faa_basename##*/}
    #    echo "  ${query_faa_basename}: ${query_faa_file}" >> ${output_config_filepath}

    #done

    # Query faa genome files (reference)
    printf "query_genome_orfs:" >> ${output_config_filepath}
    local query_ref_faa_files=($(find ${query_reference_faa_directory} -maxdepth 1 -type f -name "*.faa" | sort -h | xargs realpath))
    # TODO - confirm there is only one
    printf "${query_ref_faa_files}\n" >> ${output_config_filepath}
    #for query_ref_faa_file in ${query_ref_faa_files[@]}; do

    #    query_ref_faa_basename=${query_ref_faa_file%.faa}
    #    query_ref_faa_basename=${query_ref_faa_basename##*/}
    #    echo "  ${query_ref_faa_basename}: ${query_ref_faa_file}" >> ${output_config_filepath}

    #done

    # Subject faa genome files
    echo "# Subject files for BLAST (the name you specify will be plotted)." >> ${output_config_filepath}
    echo "  # These should be ORF predictions from the genomes of the organisms of interest." >> ${output_config_filepath}
    echo "  # If you want to plot the query genome as well, then include it here." >> ${output_config_filepath}
    echo "subjects:" >> ${output_config_filepath}
    local subject_faa_files=($(find ${subject_faa_directory} -maxdepth 1 -type f -name "*.faa" | sort -h | xargs realpath))

    for subject_faa_file in ${subject_faa_files[@]}; do

        subject_faa_basename=${subject_faa_file%.faa}
        subject_faa_basename=${subject_faa_basename##*/}
        echo "  ${subject_faa_basename}: ${subject_faa_file}" >> ${output_config_filepath}

    done

}


# Determine the run mode based on user input
# Different variables are needed for the different run modes
if [ ${run_mode} = "setup" ]; then

    # Inputs specified by the user
    # TODO - check for these
    output_config_filepath=$2
    query_faa_directory=$3
    query_reference_faa_directory=$4
    subject_faa_directory=$5

    echo "[ $(date -u) ]: Running in setup mode"

    # Generate config file
    template_config_filepath="template_config.yaml" # HARD-CODED
    cp ${template_config_filepath} ${output_config_filepath}

    add_samples_to_config_file ${output_config_filepath} ${query_faa_directory} ${query_reference_faa_directory} ${subject_faa_directory}

    # TODO - also generate templates for the user to fill in genome and gene metadata

    echo "[ $(date -u) ]: Setup complete. After modifying the config.yaml file to your liking, run BackBLAST using 'run' mode."

elif [ ${run_mode} = "run" ]; then

    # HARD-CODED variables
    snakefile=Snakefile

    # User-provided variables
    config_file=$2 # config.yaml
    run_directory=$3 # .
    jobs=$4 # 2

    snakemake --snakefile ${snakefile} --configfile ${config_file} --directory ${run_directory} --jobs ${jobs} --use-conda --reason

elif [ ${run_mode} = "auto" ]; then

    # Inputs specified by the user
    # TODO - check for these
    query_faa_directory=$2
    query_reference_faa_directory=$3
    subject_faa_directory=$4
    run_directory=$5 # .
    phylogenetic_tree_newick=$6
    jobs=$7 # 2

    # Generate config file
    template_config_filepath="template_config.yaml" # HARD-CODED
    output_config_filepath="${run_directory}/config.yaml" # HARD-CODED
    cp ${template_config_filepath} ${output_config_filepath}

    add_samples_to_config_file ${output_config_filepath} ${query_faa_directory} ${query_reference_faa_directory} ${subject_faa_directory}

    # Add phylogenetic tree
    sed -i "s/phylogenetic_tree_newick: NA/phylogenetic_tree_newick: ${phylogenetic_tree_newick}/" ${output_config_filepath}

    # HARD-CODED variables
    snakefile=Snakefile

    snakemake --snakefile ${snakefile} --configfile ${output_config_filepath} --directory ${run_directory} --jobs ${jobs}

else
    echo "[ $(date -u) ]: Provided run_mode must be 'setup', 'run', or 'auto'. You provided '${run_mode}'. Exiting..."
    exit 1
fi

echo "[ $(date -u) ]: BackBLAST finished."

