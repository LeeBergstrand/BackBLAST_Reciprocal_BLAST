#!/usr/bin/env bash
set -euo pipefail
# generate_run_templates.sh
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019
# Script for generating and manipulating BackBLAST config templates
# Part of the BackBLAST pipeline

# GLOBAL variables
# Assign ones that overlap with the main BackBLAST script only if this script is called independently
if [[ ${BASH_SOURCE[0]} = ${0} ]]; then
  readonly SCRIPT_NAME="${0##*/}"
  readonly SCRIPT_DIR="$(realpath ${0%/*})"
  readonly TEMPLATE_CONFIG="${SCRIPT_DIR}/../snakemake/template_config.yaml"
fi

#######################################
# Add the user-provided subject files to the end of the template config file
# Globals: (none)
# Arguments:
#   output_config_filepath: the path to the config.yaml file to which the subject names will be added
#   subject_genome_directory: the path to the directory containing the protein files for subject genomes (FastA format)
#   genome_extension: the extension of the predicted protein files for subject genomes (e.g., faa)
# Returns:
#   writes to output_config_filepath
#######################################
function add_subjects_to_config_file() {
  # User-provided inputs
  local output_config_filepath
  output_config_filepath=$1
  local subject_genome_directory
  subject_genome_directory=$2
  local genome_extension
  genome_extension=$3

  # Find the subjects
  local subject_genome_files
  subject_genome_files=($(find ${subject_genome_directory} -maxdepth 1 -type f -name "*.${genome_extension}" | sort -h | xargs realpath))

  echo "[ $(date -u) ]: Found ${#subject_genome_files[@]} subject genomes with extension '${genome_extension}'"

  # Append to the bottom of the file (a bit hacky)
  # Note that the arbitrary sample name is defined here for each sample as the basename of the file
  for subject_genome_file in ${subject_genome_files[@]}; do

    subject_genome_basename=${subject_genome_file%.${genome_extension}}
    subject_genome_basename=${subject_genome_basename##*/}
    echo "  ${subject_genome_basename}: '${subject_genome_file}'" >> ${output_config_filepath}

  done
}

#######################################
# Generate the genome and gene metadata file templates
# Globals: (none)
# Arguments:
#   genome_metadata_tsv: the path to which the genome metadata TSV file is to be written
#   gene_metadata_tsv: the path to which the gene metadata TSV file is to be written
#   subject_genome_directory: the path to the directory containing the protein files for subject genomes (FastA format)
#   genome_extension: the extension of the predicted protein files for subject genomes (e.g., faa)
#   query_filepath: the path to the query protein sequences (FastA format)
# Returns:
#   writes to genome_metadata_tsv and gene_metadata_tsv
#######################################
function generate_metadata_templates() {
  # Assign variables from input
  local genome_metadata_tsv
  genome_metadata_tsv=$1
  local gene_metadata_tsv
  gene_metadata_tsv=$2
  local subject_genome_directory
  subject_genome_directory=$3
  local genome_extension
  genome_extension=$4
  local query_filepath
  query_filepath=$5

  # Generate genome metadata template
  echo "[ $(date -u) ]: Writing genome metadata template to '${genome_metadata_tsv}'" >&2
  printf "subject_name\tplotting_name\n" > ${genome_metadata_tsv}

  local subject_genome_files
  subject_genome_files=($(find ${subject_genome_directory} -maxdepth 1 -type f -name "*.${genome_extension}" | sort -h))

  for subject_genome_file in ${subject_genome_files[@]}; do

    subject_genome_basename=${subject_genome_file%.${genome_extension}}
    subject_genome_basename=${subject_genome_basename##*/}
    printf "${subject_genome_basename}\t\n" >> ${genome_metadata_tsv}

  done

  # Generate gene metadata template
  echo "[ $(date -u) ]: Writing gene metadata template to '${gene_metadata_tsv}'" >&2
  printf "qseqid\tgene_name\n" > ${gene_metadata_tsv}

  local query_accessions
  query_accessions=($(grep "^>" ${query_filepath} | cut -d ">" -f 2 | cut -d " " -f 1))

  for query_accession in ${query_accessions[@]}; do

    printf "${query_accession}\t\n" >> ${gene_metadata_tsv}

  done
}

#######################################
# Generate template files for the BackBLAST run 
# Globals: (none)
# Arguments:
#   template_config: the path to the template Snakemake configuration file (YAML format)
#   query_filepath: the path to the query protein sequences (FastA format)
#   query_genome_filepath: the path to the predicted protein sequences of the entire genome corresponding to the query protein sequences (FastA format)
#   subject_genome_directory: the path to the directory containing the protein files for subject genomes (FastA format)
#   genome_extension: the extension of the predicted protein files for subject genomes (e.g., faa)
#   output_directory: path to the directory where output files should be written
#   threads: maximum number of threads that any given task within the snakemake pipeline ought to use
#   phylogenetic_tree_newick: path to the phylogenetic tree file corresponding to the subject genomes (or type 'subjects' to have the tree auto-generated)
#   bootstrap_cutoff: numeric value; only bootstrap numbers above this value (e.g., 80) will be shown on the plot (or 'NA' to skip)
#   root_name: string; the exact label of the phylogenetic tree tip corresponding to the root (or 'NA' to skip)
#   evalue: numeric/scientific; e-value cutoff for reciprocal BLASTP
#   pident: numeric; percent identity cutoff for reciprocal BLASTP
#   qcov: numeric; percent query coverage cutoff for reciprocal BLASTP
# Returns:
#   writes three template files to the output_directory: 'config.yaml', 'genome_metadata.tsv', and 'gene_metadata.tsv'
#######################################
function make_run_templates() {
  # Assign input variables
  # TODO - is there a more elegant way of doing this? This is a ton of input variables to have to bring in one at a time
  local template_config
  template_config=$1
  local query_filepath
  query_filepath=$2
  local query_genome_filepath
  query_genome_filepath=$3
  local subject_genome_directory
  subject_genome_directory=$4
  local genome_extension
  genome_extension=$5
  local output_directory
  output_directory=$6
  local threads
  threads=$7
  local phylogenetic_tree_newick
  phylogenetic_tree_newick=$8
  local bootstrap_cutoff
  bootstrap_cutoff=$9
  local root_name
  root_name=${10}
  local evalue
  evalue=${11}
  local pident
  pident=${12}
  local qcov
  qcov=${13}

  # Check if output directory exists
  if [[ ! -d ${output_directory} ]]; then
    echo "[ $(date -u) ]: ERROR: output directory '${output_directory}' does not exist. Exiting..."
    exit 1
  fi

  # Check if desired output files already exist
  local output_config_filepath
  output_config_filepath=${output_directory}/config.yaml
  local genome_metadata_tsv
  genome_metadata_tsv=${output_directory}/genome_metadata.tsv
  local gene_metadata_tsv
  gene_metadata_tsv=${output_directory}/gene_metadata.tsv

  if [[ -f ${output_config_filepath} ]]; then
    echo "[ $(date -u) ]: Found existing 'config.yaml' at '${output_directory}'. Will not continue. Exiting..." >&2
    exit 1
  elif [[ -f ${genome_metadata_tsv} ]]; then
    echo "[ $(date -u) ]: Found existing 'genome_metadata.tsv' at '${output_directory}'. Will not continue. Exiting..." >&2
    exit 1
  elif [[ -f ${gene_metadata_tsv} ]]; then
    echo "[ $(date -u) ]: Found existing 'gene_metadata.tsv' at '${output_directory}'. Will not continue. Exiting..." >&2
    exit 1
  fi

  # Generate the metadata templates
  generate_metadata_templates ${genome_metadata_tsv} ${gene_metadata_tsv} ${subject_genome_directory} \
    ${genome_extension} ${query_filepath}

  ### Generate config file and add variables
  echo "[ $(date -u) ]: Writing config info to '${output_config_filepath}'" >&2
  cp ${template_config} ${output_config_filepath}

  # Add subject info
  add_subjects_to_config_file ${output_config_filepath} ${subject_genome_directory} ${genome_extension}

  # Special check for the phylogenetic tree - if the entry is not a file (e.g., 'subjects' or 'NA'), then do not run realpath
  if [[ -f ${phylogenetic_tree_newick} ]]; then
    phylogenetic_tree_newick=$(realpath ${phylogenetic_tree_newick})
  fi

  # Add other variables to config file
  # NOTE: I'm using the '|' symbol as the sed separator because some of the variables contain forward slashes (which is the normal separator)
  # TODO - is there a more elegant way of doing this?
  sed -i "s|^query_genes: .*|query_genes: '$(realpath ${query_filepath})'|" ${output_config_filepath}
  sed -i "s|^query_genome_orfs: .*|query_genome_orfs: '$(realpath ${query_genome_filepath})'|" ${output_config_filepath}
  sed -i "s|^threads: .*|threads: ${threads}|" ${output_config_filepath}
  sed -i "s|^phylogenetic_tree_newick: .*|phylogenetic_tree_newick: '${phylogenetic_tree_newick}'|" ${output_config_filepath}
  sed -i "s|^genome_metadata_tsv: .*|genome_metadata_tsv: '${genome_metadata_tsv}'|" ${output_config_filepath}
  sed -i "s|^gene_metadata_tsv: .*|gene_metadata_tsv: '${gene_metadata_tsv}'|" ${output_config_filepath}
  sed -i "s|^bootstrap_cutoff: .*|bootstrap_cutoff: ${bootstrap_cutoff}|" ${output_config_filepath}
  sed -i "s|^root_name: .*|root_name: '${root_name}'|" ${output_config_filepath}
  sed -i "s|^e_value_cutoff: .*|e_value_cutoff: ${evalue}|" ${output_config_filepath}
  sed -i "s|^minimum_percent_identity: .*|minimum_percent_identity: ${pident}|" ${output_config_filepath}
  sed -i "s|^minimum_query_coverage: .*|minimum_query_coverage: ${qcov}|" ${output_config_filepath}
}

function main() {
  # If no input is provided, provide help and exit
  if [[ $# -eq 0 ]]; then
    echo "No arguments provided. Please run '-h' or '--help' to see help. Exiting..." >&2
    exit 1
  elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

    # Help statement
    printf "${SCRIPT_NAME}: generate template config files for a BackBLAST 'setup'. Part of the BackBLAST suite.\n"
    printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, Neufeld Research Group, 2019\n\n"
    printf "Usage: ${SCRIPT_NAME} [OPTIONS] query_filepath query_genome_filepath subject_genome_directory output_directory\n\n"
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
    printf "   -c qcov: percent query coverage cutoff for reciprocal BLASTP [default: 50]\n"
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
  local qcov
  qcov=50
  local genome_extension
  genome_extension="faa"
  local template_config
  template_config=${TEMPLATE_CONFIG}

  # Set options (help from https://wiki.bash-hackers.org/howto/getopts_tutorial; accessed March 8th, 2019)
  OPTIND=1 # reset the OPTIND counter just in case
  while getopts ":@:t:b:r:e:p:c:x:T:" opt; do
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
      c)
        qcov=${OPTARG}
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

  echo "[ $(date -u) ]: Running ${SCRIPT_NAME}" >&2
  echo "[ $(date -u) ]: Command run: ${SCRIPT_NAME} ${original_arguments}" >&2

  make_run_templates ${template_config} ${query_filepath} ${query_genome_filepath} ${subject_genome_directory} \
    ${genome_extension} ${output_directory} ${threads} ${phylogenetic_tree_newick} ${bootstrap_cutoff} ${root_name} \
    ${evalue} ${pident} ${qcov}

  echo "[ $(date -u) ]: BackBLAST template generation finished." >&2
}

# Only run the script if it is called from the command line
if [[ ${BASH_SOURCE[0]} = ${0} ]]; then
  main $@
fi

