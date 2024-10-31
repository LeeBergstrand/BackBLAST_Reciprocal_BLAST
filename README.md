BackBLAST reciprocal BLAST workflow
==========================
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3465954.svg)](https://doi.org/10.5281/zenodo.3465954)

Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2024

# Software overview
`backblast` automates the use of NCBI BLASTP to search for genes or gene clusters within a multiple bacterial genomes. 
Non-orthologous genes are filtered out by identifying and extracting only bidirectional best BLASTP hits using a graph-based algorithm.
Furthermore, `backblast` allows users to visualize the results from bidirectional BLASTP in the form of a heatmap.

The bidirectional BLASTP-based filtering algorithm is illustrated below:

![BackBLAST Algorithm](https://private-user-images.githubusercontent.com/18713012/381830418-2b8690db-ffd5-4fe5-adc3-661c6a7515c2.gif?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MzAzNTE2MDYsIm5iZiI6MTczMDM1MTMwNiwicGF0aCI6Ii8xODcxMzAxMi8zODE4MzA0MTgtMmI4NjkwZGItZmZkNS00ZmU1LWFkYzMtNjYxYzZhNzUxNWMyLmdpZj9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNDEwMzElMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjQxMDMxVDA1MDgyNlomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTVlOGJkODI1ZjIxYmE1ZTQyMTAzMzdmMjBiYjY4YzYxZGQ2ZWJlOWQ2MDkxMDY1YTE4M2U0ZWM2ODYyZmNiNmQmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.mywfr_T158k973yrDKpOsR-uHRDg_s155-a4bm2n6Hs)

Here is an example heatmap visualization of the results:

![Example Results](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41396-020-0650-2/MediaObjects/41396_2020_650_Fig7_HTML.png)
(Example from Spasov, Tsuji, _et al._, 2020, [doi:10.1038/s41396-020-0650-2](https://doi.org/10.1038/s41396-020-0650-2))

# Dependencies
- OS: Linux operating system (e.g., Ubuntu) or MacOS (tested on Sonoma 14)
- Software: miniconda or miniforge needs to be installed (must be set to the `osx-64` channel for MacOS)
- Hardware: Workflow is pretty light on RAM, CPU, and storage space, so most modern computers should be able to handle BackBLAST without issue.
  The only exception is if you create a genome tree within the pipeline, in which case you'll need a fair amount of CPU and time to calculate large trees.

# Installation
Installation is currently semi-manual. We hope to make an automated conda install in the future.

Run the following code in your command line (e.g., Terminal) to install BackBLAST:
```bash
# Download the repo
git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd BackBLAST_Reciprocal_BLAST
git checkout develop # optionally go to a specific branch or version tag

# Create the conda env
conda env create -n backblast --file=environment.yml

# Copy the key repo contents into a conda share folder
conda activate backblast
mkdir -p ${CONDA_PREFIX}/share/backblast
cp -r * ${CONDA_PREFIX}/share/backblast

# Remove the original repo
cd ..
rm -rf BackBLAST_Reciprocal_BLAST

# Add instructions to export the BackBLAST folder to the PATH when the repo activates
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

if [[ ! -f ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh ]]; then
  echo '#!/bin/sh' > ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
fi
echo "export PATH=\${PATH}:${CONDA_PREFIX}/share/backblast:${CONDA_PREFIX}/share/backblast/scripts" \
  >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
chmod 755 ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# Re-activate the repo to apply the changes
conda deactivate
conda activate backblast
```
Now you should be good to go! Run `backblast -h` to get started.

# Usage (quick start)
See the full manual + settings for BackBLAST at the bottom of this README. Quick start instructions are below.

## Recommended workflow
```bash
# Set up the run
backblast setup query.faa query_genome.faa subject_dir output_dir
# Then edit output_dir/config.yaml
# You can also edit output_dir/gene_metadata.tsv and output_dir/genome_metadata.tsv to make the plot look better

# Start the run
backblast run output_dir/config.yaml output_dir
# All done! You can iteratively refine the plot from here as you'd like.
```

## Speedy workflow
Gets the job done without any custom settings
```bash
backblast auto [OPTIONS] query.faa query_genome.faa subject_dir output_dir
```

# Test data
Try a test run from inside the repo with:
```bash
mkdir -p testing/outputs
# Make sure backblast is added to your PATH before running the test
backblast run testing/inputs/config.yaml testing/outputs --notemp

# See if the output files looks as expected
cmp testing/outputs/blast/combine_blast_tables/blast_tables_combined.csv \
  testing/outputs_expected/blast/combine_blast_tables/blast_tables_combined.csv
cmp testing/outputs/heatmap/BackBLAST_heatmap.tsv \
  testing/outputs_expected/heatmap/BackBLAST_heatmap.tsv

# Clean up test if everything looks good
rm -r testing/outputs
```

# Citation
We hope that BackBLAST is helpful for you! If you use BackBLAST in your research, please cite the following paper, which describes the initial version of BackBLAST:

> Bergstrand LH, Cardenas E, Holert J, Van Hamme JD, Mohn WW. Delineation of Steroid-Degrading Microorganisms through Comparative Genomic Analysis. __mBio__ 7:[10.1128/mbio.00166-16](https://doi.org/10.1128/mbio.00166-16).

# Full usage instructions
Help messages (e.g., from running `backblast -h`) are pasted below.

`backblast`
```
backblast: pipeline to search for and visualize gene homologs across multiple genomes.

Please specify a run mode for the main workflow for further usage instructions:
   1. setup: for setting up pre-run configuration files
   2. run: to start a run using configuration files
   3. auto: to skip setup and run the pipeline end-to-end with default values

Or run a specific step in the workflow manually:
   search: performs reciprocal BLASTP on given inputs
   remove_duplicates: removes duplicate hits from reciprocal BLASTP results
   create_blank_results: creates a blank output BLAST table for entries with no reciprocal BLASTP hits
   combine_tables: combines a set of input BLAST tables
   generate_heatmap: generates a heatmap from input BLAST and phylogeny data

Advanced parameters (use with care):
   -U utils_dir: Path to the BackBLAST utility script directory for running specific workflow steps
                 [default: /Users/jmtsuji/mambaforge/envs/backblast/share/backblast/scripts]
```

`backblast setup`
```
backblast setup: module for setting up a BackBLAST run.

Usage: backblast setup [OPTIONS] query_filepath query_genome_filepath subject_genome_directory output_directory

Positional arguments (required):
   query_filepath: path to the query predicted protein sequences from the query genome, FastA format
   query_genome_filepath: path to the predicted proteins of the entire query genome, FastA format
   subject_genome_directory: directory containing predicted proteins of all subject genomes (FastA format).
                                 One genome per file, with extension 'faa' (or specify -x)
   output_directory: directory where config ('config.yaml') and metadata templates ('genome_metadata.tsv', 'gene_metadata.tsv') will be created

Optional arguments:
   -t phylogenetic_tree_newick: path to the pre-calculated phylogenetic tree [default: 'subjects' - auto-calculate tree]
                                Specify 'NA' to skip tree generation and plot the heatmap alone.
   -b bootstrap_cutoff: numerical value (e.g., 80) under which bootstrap values will not be displayed [default: NA]
   -r root_name: Exact name of the tree tip label at the desired root of the tree [default: NA; will skip rooting]
   -e evalue: e-value cutoff for reciprocal BLASTP [default: 1e-40]
   -p pident: percent identity cutoff for reciprocal BLASTP [default: 25]
   -c qcov: percent query coverage cutoff for reciprocal BLASTP [default: 50]
   -x genome_extension: extension for predicted protein files of subject genomes [default: faa]
   -@ threads: maximum threads to use for any process [default: 1]

Advanced parameters (use with care):
   -T template_config: Path to the template config file used in setup [default: /Users/jmtsuji/mambaforge/envs/backblast/share/backblast/snakemake/template_config.yaml]
   -U utils_dir: Path to the BackBLAST utility script directory [default: /Users/jmtsuji/mambaforge/envs/backblast/share/backblast/scripts]

Note: Currently does NOT support whitespaces in any input variables.
```

`backblast run`
```
backblast run: module for initializing a BackBLAST run.

Usage: backblast run [OPTIONS] config_filepath run_directory [SNAKEMAKE_ARGUMENTS]

Positional arguments (required):
   config_filepath: path to the config file generated using the 'setup' module
   run_directory: the directory where BackBLAST results out to be output

Optional arguments:
   -P conda_prefix: path to where the conda envs should be stored [default: '.snakemake/conda' in the run_directory]
   -j jobs: Number of processing threads available for the run [default: 1]
            **Should be no lower than the 'threads' setting in the config file**

Advanced parameters (use with care):
   -S snakefile: Path to the Snakefile used to run BackBLAST [default: /Users/jmtsuji/mambaforge/envs/backblast/share/backblast/snakemake/Snakefile]
   -C use_conda: specify either 'True' or 'False' for whether or not each job should be run in 
             its own conda env [default: True]
             If 'False' is set, then all dependencies need to be installed in the main environment
             where BackBLAST is running. Could be tricky.

Snakemake arguments:
   Any flags added at the end of the command will be passed directly to snakemake, e.g., --notemp

Note: Currently does NOT support whitespaces in any input variables.
```

`backblast auto`
```
backblast auto: module for setting up a BackBLAST run AND starting with some defaults.

Usage: backblast auto [OPTIONS] query_filepath query_genome_filepath subject_genome_directory output_directory [SNAKEMAKE_ARGUMENTS]

Positional arguments (required):
   query_filepath: path to the query predicted protein sequences from the query genome, FastA format
   query_genome_filepath: path to the predicted proteins of the entire query genome, FastA format
   subject_genome_directory: directory containing predicted proteins of all subject genomes (FastA format).
                                 One genome per file, with extension 'faa' (or specify -x)
   output_directory: directory where config ('config.yaml') and metadata templates ('genome_metadata.tsv', 'gene_metadata.tsv') will be created

Optional arguments:
   -t phylogenetic_tree_newick: path to the pre-calculated phylogenetic tree [default: 'subjects' - auto-calculate tree]
                                Specify 'NA' to skip tree generation and plot the heatmap alone.
   -b bootstrap_cutoff: numerical value (e.g., 80) under which bootstrap values will not be displayed [default: NA]
   -r root_name: Exact name of the tree tip label at the desired root of the tree [default: NA; will skip rooting]
   -e evalue: e-value cutoff for reciprocal BLASTP [default: 1e-40]
   -p pident: percent identity cutoff for reciprocal BLASTP [default: 25]
   -c qcov: percent query coverage cutoff for reciprocal BLASTP [default: 50]
   -x genome_extension: extension for predicted protein files of subject genomes [default: faa]
   -P conda_prefix: path to where the conda envs should be stored [default: '.snakemake/conda' in the run_directory]   -@ threads: maximum threads to use for any process [default: 1]
   -j jobs: Number of processing threads available for the run [default: 1]
            **Should be no lower than the 'threads' setting in the config file**

Advanced parameters (use with care):
   -T template_config: Path to the template config file used in setup [default: /Users/jmtsuji/mambaforge/envs/backblast/share/backblast/snakemake/template_config.yaml]
   -S snakefile: Path to the Snakefile used to run BackBLAST [default: /Users/jmtsuji/mambaforge/envs/backblast/share/backblast/snakemake/Snakefile]
   -C use_conda: specify either 'True' or 'False' for whether or not each job should be run in 
             its own conda env [default: True]
             If 'False' is set, then all dependencies need to be installed in the main environment
             where BackBLAST is running. Could be tricky.

Snakemake arguments:
   Any flags added at the end of the command will be passed directly to snakemake, e.g., --notemp

Note: Currently does NOT support whitespaces in any input variables.
```
