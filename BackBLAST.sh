#!/usr/bin/env bash
set -euo pipefail
# BackBLAST.sh
# The entry script for running BackBLAST via command line
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019

# GLOBAL variables
readonly VERSION="0.8.0"
readonly SCRIPT_NAME="${0##*/}"
readonly SCRIPT_DIR="$(realpath ${0%/*})"
readonly TEMPLATE_CONFIG="${SCRIPT_DIR}/template_config.yaml"
readonly SNAKEFILE="${SCRIPT_DIR}/Snakefile"

# Load modules
source "${SCRIPT_DIR}/scripts/BackBLAST_generate_run_templates.sh" # the 'make_run_templates' function is called

#######################################
# Run snakemake
# Globals: (none)
# Arguments:
#   snakefile: the path to the Snakefile used to run BackBLAST
#   config_file: path to the config.yaml file containing settings for the run
#   run_directory: path to the directory where output files should be written
#   conda_prefix: the path to where conda environments are stored
#   jobs: maximum number of parallel jobs for the snakemake scheduler to run
#   use_conda: character value of either 'True' or 'False' to specify whether each job should be run in its own conda environment.
#          If 'False' is specified, then all dependencies need to be installed on your machine (e.g., in a central conda env).
#   snakemake_arguments: everything after argument 6 are flags passed directly to snakemake. May contain spaces inbetween!
# Returns:
#   output files from the BackBLAST pipeline, in the run_directory
#######################################
function run_snakemake() {
  # Get input variables
  local snakefile
  snakefile=$1
  local config_file
  config_file=$2
  local run_directory
  run_directory=$3
  local conda_prefix
  conda_prefix=$4
  local jobs
  jobs=$5
  local use_conda
  use_conda=$6
  local snakemake_arguments
  # Get everything after argument 6
  snakemake_arguments=${@:7}

  # Make sure conda_prefix is an absolute path if not the default
  if [[ "${conda_prefix}" != ".snakemake/conda" ]]; then
    conda_prefix=$(realpath "${conda_prefix}")
  fi

  # Run snakemake
  if [[ ${use_conda} == "True" ]]; then
    echo "[ $(date -u) ]: Command: snakemake --snakefile '${snakefile}' --configfile '${config_file}' --directory '${run_directory}' --conda-prefix '${conda_prefix}' --jobs ${jobs} --rerun-incomplete --reason --printshellcmds --use-conda ${snakemake_arguments}"
    snakemake --snakefile "${snakefile}" --configfile "${config_file}" --directory "${run_directory}" --conda-prefix "${conda_prefix}" --jobs "${jobs}" --rerun-incomplete --reason --printshellcmds --use-conda ${snakemake_arguments}
  elif  [[ ${use_conda} == "False" ]]; then
    # No need for --conda-prefix here; do not use
    echo "[ $(date -u) ]: Command: snakemake --snakefile '${snakefile}' --configfile '${config_file}' --directory '${run_directory}' --jobs ${jobs} --rerun-incomplete --reason --printshellcmds ${snakemake_arguments}"
    snakemake --snakefile "${snakefile}" --configfile "${config_file}" --directory "${run_directory}" --jobs ${jobs} --rerun-incomplete --reason --printshellcmds ${snakemake_arguments}
  else
    echo "[ $(date -u) ]: the 'use_conda' variable must be either 'True' or 'False'; instead, '${use_conda}' was specified. Exiting..."
    exit 1
  fi
}

#######################################
# Perform the 'setup' command
# Globals:
#   SCRIPT_NAME: the name of this script
#   VERSION: the script version
#   TEMPLATE_CONFIG: the path to the template config file used in setup (YAML format)
#   OPTIND and OPTARG: used in command line parsing as part of 'getopts'
# Arguments:
#   all command line inputs for the 'setup' module - see help statement below
# Returns:
#   runs the 'setup' command end-to-end
#######################################
function perform_setup() {
  if [[ $# -eq 0 ]]; then
    echo "setup: No arguments provided. Please run 'setup -h' or 'setup --help' for help. Exiting..." >&2
    exit 1
  elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

    # Print help statement
    printf "${SCRIPT_NAME} setup: module for setting up a BackBLAST run.\n"
    printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n"
    printf "Version: ${VERSION}\n\n"
    printf "Usage: ${SCRIPT_NAME} setup [OPTIONS] query_filepath query_genome_filepath subject_genome_directory output_directory\n\n"
    printf "Positional arguments (required):\n"
    printf "   query_filepath: path to the query predicted protein sequences from the query genome, FastA format\n"
    printf "   query_genome_filepath: path to the predicted proteins of the entire query genome, FastA format\n"
    printf "   subject_genome_directory: directory containing predicted proteins of all subject genomes (FastA format).\n"
    printf "                                 One genome per file, with extension 'faa' (or specify -x)\n"
    printf "   output_directory: directory where config ('config.yaml') and metadata templates ('genome_metadata.tsv', 'gene_metadata.tsv') will be created\n\n"
    printf "Optional arguments:\n"
    printf "   -t phylogenetic_tree_newick: path to the pre-calculated phylogenetic tree [default: 'subjects' - auto-calculate tree]\n"
    printf "                                Specify 'NA' to skip tree generation and plot the heatmap alone.\n"
    printf "   -b bootstrap_cutoff: numerical value (e.g., 80) under which bootstrap values will not be displayed [default: NA]\n"
    printf "   -r root_name: Exact name of the tree tip label at the desired root of the tree [default: NA; will skip rooting]\n"
    printf "   -e evalue: e-value cutoff for reciprocal BLASTP [default: 1e-40]\n"
    printf "   -p pident: percent identity cutoff for reciprocal BLASTP [default: 25]\n"
    printf "   -x genome_extension: extension for predicted protein files of subject genomes [default: faa]\n"
    printf "   -@ threads: maximum threads to use for any process [default: 1]\n\n"
    printf "Advanced parameters (use with care):\n"
    printf "   -T template_config: Path to the template config file used in setup [default: ${TEMPLATE_CONFIG}]\n\n"
    printf "Note: Currently does NOT support whitespaces in any input variables.\n\n"

    # Exit
    exit 0
  fi

  # Set defaults for options
  local threads
  threads=1
  local phylogenetic_tree_newick
  phylogenetic_tree_newick="subjects"
  local bootstrap_cutoff
  bootstrap_cutoff="NA"
  local root_name
  root_name="NA"
  local evalue
  evalue=1e-40
  local pident
  pident=25
  local genome_extension
  genome_extension="faa"
  local template_config
  template_config=${TEMPLATE_CONFIG}

  # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
  OPTIND=1 # reset the OPTIND counter just in case
  while getopts ":@:t:b:r:e:p:x:T:" opt; do
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
      x)
        genome_extension=${OPTARG}
        ;;
      T)
        template_config=${OPTARG}
        ;;
      \?)
        echo "[ $(date -u) ]: ERROR: Invalid option: '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
      :)
        echo "[ $(date -u) ]: ERROR: argument needed following '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
    esac
  done

  # Set positional arguments
  local original_arguments
  original_arguments=${@} # save for reporting later
  shift $((OPTIND - 1)) # shift to avoid flags when assigning positional arguments
  local query_filepath
  query_filepath=$1
  local query_genome_filepath
  query_genome_filepath=$2
  local subject_genome_directory
  subject_genome_directory=$3
  local output_directory
  output_directory=$4

  echo "[ $(date -u) ]: Running ${SCRIPT_NAME} in setup mode" >&2
  echo "[ $(date -u) ]: Command run: ${SCRIPT_NAME} setup ${original_arguments}" >&2

  # Make the template files for the run
  make_run_templates ${template_config} ${query_filepath} ${query_genome_filepath} ${subject_genome_directory} \
    ${genome_extension} ${output_directory} ${threads} ${phylogenetic_tree_newick} ${bootstrap_cutoff} ${root_name} ${evalue} ${pident}

  echo "[ $(date -u) ]: Setup complete. After modifying the config.yaml file to your liking, run ${SCRIPT_NAME} using 'run' mode." >&2
}

#######################################
# Perform the 'run' command
# Globals:
#   SCRIPT_NAME: the name of this script
#   VERSION: the script version
#   SNAKEFILE: the path to the Snakefile used to run BackBLAST
# Arguments:
#   all command line inputs for the 'run' module - see help statement below
# Returns:
#   runs the 'run' command end-to-end
#######################################
function perform_run() {

  if [[ $# -eq 0 ]]; then
    echo "run: No arguments provided. Please run 'setup -h' or 'setup --help' for help. Exiting..." >&2
    exit 1
  elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

    # Print help statement
    printf "${SCRIPT_NAME} run: module for initializing a BackBLAST run.\n"
    printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n"
    printf "Version: ${VERSION}\n\n"
    printf "Usage: ${SCRIPT_NAME} run [OPTIONS] config_filepath run_directory [SNAKEMAKE_ARGUMENTS]\n\n"
    printf "Positional arguments (required):\n"
    printf "   config_filepath: path to the config file generated using the 'setup' module\n"
    printf "   run_directory: the directory where BackBLAST results out to be output\n\n"
    printf "Optional arguments:\n"
    printf "   -P conda_prefix: path to where the conda envs should be stored [default: '.snakemake/conda' in the run_directory]\n"
    printf "   -j jobs: Number of processing threads available for the run [default: 1]\n"
    printf "            **Should be no lower than the 'threads' setting in the config file**\n\n"
    printf "Advanced parameters (use with care):\n"
    printf "   -S snakefile: Path to the Snakefile used to run BackBLAST [default: ${SNAKEFILE}]\n"
    printf "   -C use_conda: specify either 'True' or 'False' for whether or not each job should be run in \n"
    printf "             its own conda env [default: True]\n"
    printf "             If 'False' is set, then all dependencies need to be installed in the main environment\n"
    printf "             where BackBLAST is running. Could be tricky.\n\n"
    printf "Snakemake arguments:\n"
    printf "   Any flags added at the end of the command will be passed directly to snakemake, e.g., --notemp\n\n"
    printf "Note: Currently does NOT support whitespaces in any input variables.\n\n"

    # Exit
    exit 0
  fi

  # Set defaults for options
  local jobs
  jobs=1
  local conda_prefix
  conda_prefix=".snakemake/conda"
  local snakefile
  snakefile=${SNAKEFILE}
  local use_conda
  use_conda="True"

  # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
  OPTIND=1 # reset the OPTIND counter just in case
  while getopts ":j:P:S:C:" opt; do
    case ${opt} in
      j)
        jobs=${OPTARG}
        ;;
      P)
        conda_prefix=${OPTARG}
        ;;
      S)
        snakefile=${OPTARG}
        ;;
      C)
        use_conda=${OPTARG}
        ;;
      \?)
        echo "[ $(date -u) ]: ERROR: Invalid option: '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
      :)
        echo "[ $(date -u) ]: ERROR: argument needed following '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
    esac
  done

  # Set positional arguments
  local original_arguments
  original_arguments=${@} # save for reporting later
  shift $((OPTIND - 1)) # shift to avoid flags when assigning positional arguments
  local config_filepath
  config_filepath=$1 # config.yaml
  local run_directory
  run_directory=$2

  # Get Snakemake arguments; i.e., everything after argument 2
  local snakemake_arguments
  snakemake_arguments=${@:3}

  echo "[ $(date -u) ]: Running ${SCRIPT_NAME} in run mode" >&2
  echo "[ $(date -u) ]: Command run: ${SCRIPT_NAME} run ${original_arguments}" >&2

  # Start the run
  run_snakemake ${snakefile} ${config_filepath} ${run_directory} ${conda_prefix} ${jobs} ${use_conda} ${snakemake_arguments}
}

#######################################
# Perform the 'auto' command
# Globals:
#   SCRIPT_NAME: the name of this script
#   VERSION: the script version
#   TEMPLATE_CONFIG: the path to the template config file used in setup (YAML format)
#   SNAKEFILE: the path to the Snakefile used to run BackBLAST
#   OPTIND and OPTARG: used in command line parsing as part of 'getopts'
# Arguments:
#   all command line inputs for the 'auto' module - see help statement below
# Returns:
#   runs the 'auto' command end-to-end (like running 'setup', then 'run')
#######################################
# TODO - consider allowing the user to specify pre-made genome and gene metadata files
function perform_auto() {
  if [[ $# -eq 0 ]]; then
    echo "auto: No arguments provided. Please run 'auto -h' or 'auto --help' for help. Exiting..." >&2
    exit 1
  elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

    # Print help statement
    printf "${SCRIPT_NAME} auto: module for setting up a BackBLAST run AND starting with some defaults.\n"
    printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n"
    printf "Version: ${VERSION}\n\n"
    printf "Usage: ${SCRIPT_NAME} auto [OPTIONS] query_filepath query_genome_filepath subject_genome_directory output_directory [SNAKEMAKE_ARGUMENTS]\n\n"
    printf "Positional arguments (required):\n"
    printf "   query_filepath: path to the query predicted protein sequences from the query genome, FastA format\n"
    printf "   query_genome_filepath: path to the predicted proteins of the entire query genome, FastA format\n"
    printf "   subject_genome_directory: directory containing predicted proteins of all subject genomes (FastA format).\n"
    printf "                                 One genome per file, with extension 'faa' (or specify -x)\n"
    printf "   output_directory: directory where config ('config.yaml') and metadata templates ('genome_metadata.tsv', 'gene_metadata.tsv') will be created\n\n"
    printf "Optional arguments:\n"
    printf "   -t phylogenetic_tree_newick: path to the pre-calculated phylogenetic tree [default: 'subjects' - auto-calculate tree]\n"
    printf "                                Specify 'NA' to skip tree generation and plot the heatmap alone.\n"
    printf "   -b bootstrap_cutoff: numerical value (e.g., 80) under which bootstrap values will not be displayed [default: NA]\n"
    printf "   -r root_name: Exact name of the tree tip label at the desired root of the tree [default: NA; will skip rooting]\n"
    printf "   -e evalue: e-value cutoff for reciprocal BLASTP [default: 1e-40]\n"
    printf "   -p pident: percent identity cutoff for reciprocal BLASTP [default: 25]\n"
    printf "   -x genome_extension: extension for predicted protein files of subject genomes [default: faa]\n"
    printf "   -P conda_prefix: path to where the conda envs should be stored [default: '.snakemake/conda' in the run_directory]"
    printf "   -@ threads: maximum threads to use for any process [default: 1]\n"
    printf "   -j jobs: Number of processing threads available for the run [default: 1]\n"
    printf "            **Should be no lower than the 'threads' setting in the config file**\n\n"
    printf "Advanced parameters (use with care):\n"
    printf "   -T template_config: Path to the template config file used in setup [default: ${TEMPLATE_CONFIG}]\n"
    printf "   -S snakefile: Path to the Snakefile used to run BackBLAST [default: ${SNAKEFILE}]\n"
    printf "   -C use_conda: specify either 'True' or 'False' for whether or not each job should be run in \n"
    printf "             its own conda env [default: True]\n"
    printf "             If 'False' is set, then all dependencies need to be installed in the main environment\n"
    printf "             where BackBLAST is running. Could be tricky.\n\n"
    printf "Snakemake arguments:\n"
    printf "   Any flags added at the end of the command will be passed directly to snakemake, e.g., --notemp\n\n"
    printf "Note: Currently does NOT support whitespaces in any input variables.\n\n"

    # Exit
    exit 0
  fi

  # Set defaults for options
  local jobs
  jobs=1
  local threads
  threads=1
  local conda_prefix
  conda_prefix=".snakemake/conda"
  local phylogenetic_tree_newick
  phylogenetic_tree_newick="subjects"
  local bootstrap_cutoff
  bootstrap_cutoff="NA"
  local root_name
  root_name="NA"
  local evalue
  evalue=1e-40
  local pident
  pident=25
  local genome_extension
  genome_extension="faa"
  local template_config
  template_config=${TEMPLATE_CONFIG}
  local snakefile
  snakefile=${SNAKEFILE}
  local use_conda
  use_conda="True"

  # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
  OPTIND=1 # reset the OPTIND counter just in case
  while getopts ":j:@:P:t:b:r:e:p:x:T:S:C:" opt; do
    case ${opt} in
      j)
        jobs=${OPTARG}
        ;;
      \@)
        threads=${OPTARG}
        ;;
      P)
        conda_prefix=${OPTARG}
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
      x)
        genome_extension=${OPTARG}
        ;;
      T)
        template_config=${OPTARG}
        ;;
      S)
        snakefile=${OPTARG}
        ;;
      C)
        use_conda=${OPTARG}
        ;;
      \?)
        echo "[ $(date -u) ]: ERROR: Invalid option: '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
      :)
        echo "[ $(date -u) ]: ERROR: argument needed following '-${OPTARG}'. Exiting..." >&2
        exit 1
        ;;
    esac
  done

  # Set positional arguments
  local original_arguments
  original_arguments=${@} # save for reporting later
  shift $((OPTIND - 1)) # shift to avoid flags when assigning positional arguments
  local query_filepath
  query_filepath=$1
  local query_genome_filepath
  query_genome_filepath=$2
  local subject_genome_directory
  subject_genome_directory=$3
  local output_directory
  output_directory=$4

  # Get Snakemake arguments; i.e., everything after argument 4
  local snakemake_arguments
  snakemake_arguments=${@:5}

  echo "[ $(date -u) ]: Running ${SCRIPT_NAME} in 'auto' mode" >&2
  echo "[ $(date -u) ]: Command run: ${SCRIPT_NAME} setup ${original_arguments}" >&2

  # Make the template files for the run
  make_run_templates ${template_config} ${query_filepath} ${query_genome_filepath} ${subject_genome_directory} \
    ${genome_extension} ${output_directory} ${threads} ${phylogenetic_tree_newick} ${bootstrap_cutoff} ${root_name} ${evalue} ${pident}

  # Change metadata files to NA for default run
  local output_config_filepath
  output_config_filepath=${output_directory}/config.yaml
  sed -i "s|^genome_metadata_tsv: .*|genome_metadata_tsv: NA|" ${output_config_filepath}
  sed -i "s|^gene_metadata_tsv: .*|gene_metadata_tsv: NA|" ${output_config_filepath}

  # Start the run
  echo "[ $(date -u) ]: Default setup finished; starting pipeline" >&2
  run_snakemake ${snakefile} ${output_config_filepath} ${output_directory} ${conda_prefix} ${jobs} ${use_conda} ${snakemake_arguments}
}

function main() {
  # If no input is provided, provide help and exit
  if [[ $# -eq 0 ]]; then
    echo "No arguments provided. Please run '-h' or '--help' to see help. Exiting..." >&2
    exit 1
  elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

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

  # Get user-supplied variables
  local run_mode
  run_mode=$1

  # Determine the run mode based on user input
  # Different variables are needed for the different run modes
  # Only pass the second and upward arguments via ${@:2} - see https://stackoverflow.com/a/9057392 (accessed July 31, 2019)
  if [[ ${run_mode} = "setup" ]]; then
    perform_setup ${@:2}
  elif [[ ${run_mode} = "run" ]]; then
    perform_run ${@:2}
  elif [[ ${run_mode} = "auto" ]]; then
    perform_auto ${@:2}
  else
    echo "Provided run_mode must be 'setup', 'run', or 'auto'. You provided '${run_mode}'. Exiting..." >&2
    exit 1
  fi

  echo "[ $(date -u) ]: BackBLAST finished." >&2
}

# Only run the script if it is called from the command line
if [[ ${BASH_SOURCE[0]} = ${0} ]]; then
  main $@
fi

