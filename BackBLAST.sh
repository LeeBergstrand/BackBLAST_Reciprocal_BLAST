#!/usr/bin/env bash
set -euo pipefail
# BackBLAST.sh
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019
# This is the entry script for running BackBLAST

# HARD-CODED variables
VERSION=1.9.9-alpha
SCRIPT_NAME=${0##*/}
SCRIPT_NAME=${SCRIPT_NAME%.*}
TEMPLATE_CONFIG="template_config.yaml"
SNAKEFILE=Snakefile

# If no input is provided, provide help and exit
if [ $# -eq 0 ]; then

    echo "No arguments provided. Please run '-h' or '--help' to see help. Exiting..."
    exit 1

elif [ $1 = "-h" -o $1 = "--help" ]; then

    # Help statement
    printf "${SCRIPT_NAME}: pipeline to search for and visualize gene homologs across multiple genomes.\n"
    printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n"
    printf "Version: ${VERSION}\n\n"
    printf "Please specify a run mode for further usage instructions:\n"
    printf "   1. setup: for setting up pre-run configuration files\n"
    printf "   2. run: to start a run using configuration files\n"
    printf "   3. auto: to skip setup and run the pipeline end-to-end with default values\n\n"

    # Exit
    exit 0

fi


# Function to add the user-provided samples to the end of the template config file
function add_samples_to_config_file {
    # User-provided inputs
    output_config_filepath=$1
    query_faa_filepath=$2
    query_faa_genome_filepath=$3
    subject_faa_directory=$4

    # Query files
    # TODO - consider moving these to main setup function
    echo "query_genes: $(realpath ${query_faa_filepath})" >> ${output_config_filepath}
    echo "query_genome_orfs: $(realpath ${query_faa_genome_filepath})" >> ${output_config_filepath}

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

# Function to run snakemake
function run_snakemake {

    # Get input variables
    config_file=$1
    run_directory=$2
    jobs=$3

    snakemake --snakefile ${SNAKEFILE} --configfile ${config_file} --directory ${run_directory} --jobs ${jobs} --use-conda --reason

}

# Function to perform the setup command
function perform_setup {

    if [ $# -eq 0 ]; then
        echo "setup: No arguments provided. Please run 'setup -h' or 'setup --help' for help. Exiting..."
        exit 1
    elif [ $1 = "-h" -o $1 = "--help" ]; then

        # Print help statement
        printf "${SCRIPT_NAME} setup: module for setting up a BackBLAST run.\n"
        printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n"
        printf "Version: ${VERSION}\n\n"
        printf "Usage: ${0##*/} setup [OPTIONS] query_faa_filepath query_faa_genome_filepath subject_faa_genome_directory output_directory\n\n"
        printf "Positional arguments (required):\n"
        printf "   query_faa_filepath: path to the query predicted protein sequences from the query genome, with extension '.faa'\n"
        printf "   query_faa_genome_filepath: path to the predicted proteins of the entire query genome, with extension '.faa'\n"
        printf "   subject_faa_genome_directory: directory containing predicted proteins of all subject genomes. One genome per file, with extension '.faa'\n"
        printf "   output_directory: directory where config ('config.yaml') and metadata templates ('genome_metadata.tsv', 'gene_metadata.tsv') will be created\n\n"
        printf "Optional arguments:\n"
        printf "   -t phylogenetic_tree_newick: path to the pre-calculated phylogenetic tree [default: 'subjects' - auto-calculate tree]\n"
        printf "   -b bootstrap_cutoff: numerical value (e.g., 80) under which bootstrap values will not be displayed [default: NA]\n"
        printf "   -r root_name: Exact name of the tree tip label at the desired root of the tree [default: NA; will midpoint root]\n"
        printf "   -e evalue: e-value cutoff for reciprocal BLASTP [default: 1e-40]\n"
        printf "   -p pident: percent identity cutoff for reciprocal BLASTP [default: 25]\n"
        printf "   -@ threads: maximum threads to use for any process [default: 1]\n\n"

        # Exit
        exit 0

    fi

    # Set defaults for options
    local threads=1
    local phylogenetic_tree_newick="subjects"
    local bootstrap_cutoff="NA"
    local root_name="NA"
    local evalue=1e-40
    local pident=25

    # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
    OPTIND=1 # reset the OPTIND counter just in case
    while getopts ":@:t:b:r:e:p:" opt; do
	    case ${opt} in
            \@)
                threads=${OPTARG}
                ;;
            t)
                phylogenetic_tree_newick=${OPTARG}
                ;;
            b)
                bootstrap_cutoff=${OPTARG}
                ;;
            r)
                root_name=${OPTARG}
                ;;
            e)
                evalue=${OPTARG}
                ;;
            p)
                pident=${OPTARG}
                ;;
		    \?)
			    (>&2 echo "[ $(date -u) ]: ERROR: Invalid option: '-${OPTARG}'. Exiting...")
			    exit 1
			    ;;
		    :)
			    (>&2 echo "[ $(date -u) ]: ERROR: argument needed following '-${OPTARG}'. Exiting...")
			    exit 1
			    ;;
	    esac
    done

    # Set positional arguments
    local original_arguments=${@} # save for reporting later
    shift $((OPTIND - 1)) # shift to avoid flags when assigning positional arguments
    local query_faa_filepath=$1
    local query_faa_genome_filepath=$2
    local subject_faa_directory=$3
    local output_directory=$4

    echo "[ $(date -u) ]: Running BackBLAST in setup mode"
    echo "[ $(date -u) ]: Command run: ${0##*/} setup ${original_arguments}"

    # Check if desired output files already exist
    local output_config_filepath=${output_directory}/config.yaml
    local genome_metadata_tsv=${output_directory}/genome_metadata.tsv
    local gene_metadata_tsv=${output_directory}/gene_metadata.tsv

    if [ -f ${output_config_filepath} ]; then
        echo "[ $(date -u) ]: Found existing 'config.yaml' at '${output_directory}'. Will not continue. Exiting..."
        exit 1
    elif [ -f ${genome_metadata_tsv} ]; then
        echo "[ $(date -u) ]: Found existing 'genome_metadata.tsv' at '${output_directory}'. Will not continue. Exiting..."
        exit 1
    elif [ -f ${gene_metadata_tsv} ]; then
        echo "[ $(date -u) ]: Found existing 'gene_metadata.tsv' at '${output_directory}'. Will not continue. Exiting..."
        exit 1
    fi

    # Generate config file and add samples
    cp ${TEMPLATE_CONFIG} ${output_config_filepath}
    add_samples_to_config_file ${output_config_filepath} ${query_faa_filepath} ${query_faa_genome_filepath} ${subject_faa_directory}

    # Add other variables
    # NOTE: I'm using the '|' symbol as the sed separator because some of the variables contain forward slashes (which is the normal separator)
    # TODO - is there a more elegant way of doing this?
    sed -i "s|^threads: .*|threads: ${threads}|" ${output_config_filepath}
    sed -i "s|^phylogenetic_tree_newick: .*|phylogenetic_tree_newick: $(realpath ${phylogenetic_tree_newick})|" ${output_config_filepath}
    sed -i "s|^genome_metadata_tsv: .*|genome_metadata_tsv: $(realpath ${genome_metadata_tsv})|" ${output_config_filepath}
    sed -i "s|^gene_metadata_tsv: .*|gene_metadata_tsv: $(realpath ${gene_metadata_tsv})|" ${output_config_filepath}
    sed -i "s|^bootstrap_cutoff: .*|bootstrap_cutoff: ${bootstrap_cutoff}|" ${output_config_filepath}
    sed -i "s|^root_name: .*|root_name: ${root_name}|" ${output_config_filepath}
    sed -i "s|^e_value_cutoff: .*|e_value_cutoff: ${evalue}|" ${output_config_filepath}
    sed -i "s|^minimum_percent_identity: .*|minimum_percent_identity: ${pident}|" ${output_config_filepath}

    # Generate genome metadata template
    echo "[ $(date -u) ]: Writing genome metadata template to '${genome_metadata_tsv}'"
    printf "subject_name\tplotting_name\n" > ${genome_metadata_tsv}

    local subject_faa_files=($(find ${subject_faa_directory} -maxdepth 1 -type f -name "*.faa" | sort -h))

    for subject_faa_file in ${subject_faa_files[@]}; do

        subject_faa_basename=${subject_faa_file%.faa}
        subject_faa_basename=${subject_faa_basename##*/}
        printf "${subject_faa_basename}\t\n" >> ${genome_metadata_tsv}

    done

    # Generate gene metadata template
    echo "[ $(date -u) ]: Writing gene metadata template to '${gene_metadata_tsv}'"
    printf "qseqid\tgene_name\n" > ${gene_metadata_tsv}

    local query_accessions=($(grep "^>" ${query_faa_filepath} | cut -d ">" -f 2 | cut -d " " -f 1))

    for query_accession in ${query_accessions[@]}; do

        printf "${query_accession}\t\n" >> ${gene_metadata_tsv}

    done

}

# Function to perform the 'run' command
function perform_run {

    if [ $# -eq 0 ]; then
        echo "run: No arguments provided. Please run 'setup -h' or 'setup --help' for help. Exiting..."
        exit 1
    elif [ $1 = "-h" -o $1 = "--help" ]; then

        # Print help statement
        printf "${SCRIPT_NAME} run: module for initializing a BackBLAST run.\n"
        printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n"
        printf "Version: ${VERSION}\n\n"
        printf "Usage: ${0##*/} run config_filepath run_directory jobs\n\n"
        printf "Positional arguments (required):\n"
        printf "   config_filepath: path to the config file generated using the 'setup' module\n"
        printf "   run_directory: the directory where BackBLAST results out to be output\n"
        printf "   jobs: Number of processing threads available for the run\n\n"

        # Exit
        exit 0

    fi

    # User-provided variables
    config_file=$1 # config.yaml
    run_directory=$2
    jobs=$3

    # Start the run
    run_snakemake ${config_file} ${run_directory} ${jobs}

}

# Function to perform the 'auto' command
# TODO - this needs to be refactored for clarity. Initial working version only.
function perform_auto {

    # TODO - make this more informative and figure out how to separate from 'setup' parser
    if [ $# -eq 0 ]; then
        echo "auto: No arguments provided. Please run 'auto -h' or 'auto --help' for help. Exiting..."
        exit 1
    elif [ $1 = "-h" -o $1 = "--help" ]; then

        # Print help statement
        printf "############################################\n"
        printf "${SCRIPT_NAME} auto: module for setting up a BackBLAST run AND starting with some defaults.\n"
        printf "Usage: same as 'setup' module, except that a run immediately proceeds afterwards. Setup module help is below.\n"
        printf "############################################\n\n"

        # DO NOT EXIT - continue to setup module

    fi

    # First set up the run
    # TODO - make a dedicated help statement for this function
    perform_setup ${@}

    # Then change metadata files to NA
    local output_directory=$4
    local output_config_filepath=${output_directory}/config.yaml
    sed -i "s|^genome_metadata_tsv: .*|genome_metadata_tsv: NA|" ${output_config_filepath}
    sed -i "s|^gene_metadata_tsv: .*|gene_metadata_tsv: NA|" ${output_config_filepath}

    # Get jobs
    # TODO - clarify that 'threads' = 'jobs'
    local jobs=$(grep "^threads: " ${output_config_filepath} | cut -d " " -f 2) # TODO - a bit hacky

    # Start the run
    run_snakemake ${config_file} ${run_directory} ${jobs}

}

function main {
    run_mode=$1

    # Determine the run mode based on user input
    # Different variables are needed for the different run modes
    if [ ${run_mode} = "setup" ]; then
        perform_setup ${@:2}
        echo "[ $(date -u) ]: Setup complete. After modifying the config.yaml file to your liking, run BackBLAST using 'run' mode."
    elif [ ${run_mode} = "run" ]; then
        perform_run ${@:2}
    elif [ ${run_mode} = "auto" ]; then
        perform_auto ${@:2}
    else
        echo "Provided run_mode must be 'setup', 'run', or 'auto'. You provided '${run_mode}'. Exiting..."
        exit 1
    fi

    echo "[ $(date -u) ]: BackBLAST finished."
}

main $@

