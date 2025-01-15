BackBLAST reciprocal BLAST workflow
==========================
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3465954.svg)](https://doi.org/10.5281/zenodo.3465954)

Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2025

# Software overview
`backblast` automates the use of NCBI BLASTP to search for genes or gene clusters within bacterial genomes. 
Non-orthologous genes are filtered out by identifying and extracting only bidirectional best BLASTP hits using a graph-based algorithm.
`backblast` then visualizes the results of bidirectional BLASTP in a convenient gene heatmap coupled to a genome phylogeny.

The bidirectional BLASTP-based filtering algorithm is illustrated below:

![BackBLAST Algorithm](https://private-user-images.githubusercontent.com/18713012/381830418-2b8690db-ffd5-4fe5-adc3-661c6a7515c2.gif?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MzY4MzM2MTUsIm5iZiI6MTczNjgzMzMxNSwicGF0aCI6Ii8xODcxMzAxMi8zODE4MzA0MTgtMmI4NjkwZGItZmZkNS00ZmU1LWFkYzMtNjYxYzZhNzUxNWMyLmdpZj9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTAxMTQlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwMTE0VDA1NDE1NVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTRmOWM2ODMyMjlmMjRjNDRiYzI0MzRkZGVjZDZmMTYzNmE2OWYwYjRjYjg2ZGZhMTdjYTA5ZjU4NjRjN2FkNzcmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0._Zn0HWRNK6News7h6oTw2r7_GnXWOjGk4wsuqBH6fS4)

Example gene heatmap visualization:

![Example Results](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41396-020-0650-2/MediaObjects/41396_2020_650_Fig7_HTML.png)
(From Spasov, Tsuji, _et al._, 2020, [doi:10.1038/s41396-020-0650-2](https://doi.org/10.1038/s41396-020-0650-2))

# Requirements and dependencies
- OS: runs on linux (e.g., Ubuntu) and MacOS (tested on Sonoma 14 and Sequoia 15)
- Hardware: most modern computers (e.g., with basic CPU, >=4 GB RAM, and >=4 GB of free storage) should be able to run BackBLAST without issue.
  The only exception is if you create a genome tree within the pipeline, in which case you'll need a fair amount of CPU and time to calculate large trees.
- Software: miniconda or miniforge needs to be installed (must be set to the `osx-64` channel for MacOS).
  (Because conda is used during the run to install additional tools for the workflow, you will need an internet connection the first time you run BackBLAST.)

# Installation
Installation is currently semi-manual. We hope to make an automated conda install in the future.

Run the following code in your command line (e.g., Terminal) to install BackBLAST:
```bash
# Download the repo
git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd BackBLAST_Reciprocal_BLAST
# git checkout develop # optionally go to a specific branch
# git checkout v2.0.0-beta1 # optionally go to a specific release

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

# Usage: quick start
```bash
# Set up the run
backblast setup query.faa query_genome.faa subject_dir output_dir
# Then edit output_dir/config.yaml
# You can also edit output_dir/gene_metadata.tsv and output_dir/genome_metadata.tsv to make the plot look better

# Start the run
backblast run output_dir/config.yaml output_dir
# All done! You can iteratively refine the plot from here as you'd like.
```
See more detailed instructions below.

# Usage: step-by-step instructions
Get beautiful heatmaps using BackBLAST by following the step-by-step instructions below.

## 1. Prepare your input files
Prepare the following files and folders as input for BackBLAST:
- `query.faa`: multi-fasta file containing the query proteins that you want to search for. 
               These should be copied and pasted from the `query_genome.faa` file below (the FastA headers for each protein need to be the same as in that file).
- `query_genome.faa`: multi-fasta file containing all proteins encoded by the query bacterial genome.
- `subject_dir`: a folder containing the subjects for the search.
                 Each subject should be the protein-coding genes (as amino acids) for a bacterial genome.
                 Save each subject as multi-fasta file (like `query_genome.faa` above).
                 By default, the extension of each file needs to be `.faa`, but this can be changed in the tool settings if desired.
- `output_dir`: an empty folder for saving the output of BackBLAST.

You are then ready to run BackBLAST, using either method 2a or 2b below.

## 2a. Recommended workflow
Gives you the ability to customize the parameters of your run.

### 2a-i. Set up the run
```bash
backblast setup query.faa query_genome.faa subject_dir output_dir
```
This makes several configuration files in the `output_dir`.

### 2a-ii. Configure the run settings and metadata
Open and edit `output_dir/config.yaml` (e.g., in a text editor) to customize the settings for your run. Some key settings include:
- Thresholds for bidirectional BLASTP: you can set the e-value, percent identity, and query coverage cutoffs for the BLASTP search.
- Phylogenetic tree: you can also set whether you want to auto-generate a genome-based phylogenetic tree during your run (takes time)
  or if you want to supply a custom pre-generated phylogenetic tree.
  - If you provide your own phylogenetic tree, the genome IDs in the tree need to match the file names of the genomes in the `subject_dir`.
    In addition, all genomes in the phylogenetic tree need to be present in the `subject_dir`, or else BackBLAST will fail. If some
    genomes are in the `subject_dir` but are missing in the tree, then BackBLAST will still proceed but will drop those genomes from the
    analysis.)
  - Alternatively, you can also choose to not add a phylogenetic tree to the heatmap.
- Tree rooting: one other helpful setting is related to tree rooting. You can choose to root the tree automatically by midpoint or provide
  the name of one of the genomes that you want to serve as the root of the tree. Alternatively, you can choose to use the current tree
  topology as-is.
- Plot width and height: you can adjust the final width and height of the heatmap as well. It might take a few rounds of running BackBLAST
  to find the "perfect" width and height for your heatmap (see below).
- Take a look at the full config file to see other advanced settings.

You also have the opportunity (optionally) to edit the tab-separeted tables `output_dir/gene_metadata.tsv` and `output_dir/genome_metadata.tsv`
to improve data visualization. You can open these in a text editor or a table editor (like Excel).
- In `gene_metadata.tsv`, you can provide human-readable names for the query genes you want to search for. You can also customize the order
  that the genes will be shown in the heatmap by changing how they are sorted in this file. If you delete rows, then those genes will be removed
  removed from the heatmap during the run.
- In `genome_metadata.tsv`, you can provide human-readable names for the genomes in the run. The sort order of the genomes in this file only
  impacts the order that the genomes are plotted in if you choose to not use a phylogenetic tree (in the config file above). Otherwise, the
  sort order is determined using the phylogenetic tree.

### 2a-iii. Start the run
```bash
backblast run output_dir/config.yaml output_dir
```
This should generate a series of blast and heatmap-related files in the output folder.

### 2a-iv. Check and iteratively improve the run results
Then, take a look at the key output files `heatmap/BackBLAST_heatmap.tsv` and `heatmap/BackBLAST_heatmap.pdf`.
The first of these files shows the final bidirectional BLASTP results as a tab-separated table. The second is the final heatmap.

If you aren't satisfied with the results, you can selectively delete run files from the output folder, tweak the settings files above, and
then re-run BackBLAST (via step 2a-iii) to re-generate the desired files with the new run settings.

Three common things to change are:
- If your BLASTP results are too stringent or not stringent enough: delete the `blast` folder, tweak the key BLASTP settings (in the config file),
  and then re-run `backblast run`. This will redo everything from the `blast` step and will overwrite the old heatmap files.
- If your visualization doesn't look right (e.g., the height and width aren't optimal, the genome/gene names aren't ideal, or the genes aren't sorted
  as you'd like), then delete the `heatmap` folder, tweak the config and metadata files (mentioned above), and re-run `backblast run`.
  BackBLAST will then use the old `blast` results and just generate a new heatmap.
- If you want to use a different phylogenetic tree, then point to a new phylogenetic tree file in the config file, then re-run `backblast run`.
  BackBLAST will use the old `blast` results but should overwrite the old heatmap files with new ones using the updated tree.

### 2a-v. Final tweaks to the heatmap
The PDF heatmap will not be perfect (e.g., the dashed lines between the genome names and tree tips won't be perfectly aligned), but
you can open the heatmap in a PDF editor like Inkscape and clean it up or customize it as desired.

## 2b. Alternative, speedy workflow
Gets the job done without any custom settings.
This command skips the setup setup (above) and just runs the pipeline. You can set some of the key config settings via optional flags. 
```bash
backblast auto [OPTIONS] query.faa query_genome.faa subject_dir output_dir
# See the different OPTIONS for customizing your run in the full settings at the bottom of this README.
```

The level of fine-tuning possible using this method is much less than in Method 2a, but it is a nice way to quickly see some initial
bidirectional BLASTP results.

In addition, you can iteratively tweak the run results and re-run using `backblast run` just like described in step 2a-iv above.

Once done, you can then edit the PDF heatmap as described in step 2a-v above.

# Caveats
One disadvantage of BackBLAST is that it has trouble handling paralogous genes (e.g., gene duplicates).
- If the __query genome__ contains paralogs of your query gene, this can sometimes cause BackBLAST to miss real functional gene hits
  for your query gene in some subject genomes. If the gene in the subject genome happens to be a closer match to one of the paralogs
  in the query genome than to your query, then even if the gene in the subject genome is a valid functional gene for the process you are
  interested in, it will be screened out by the BackBLAST algorithm and not reported. This can lead to an under-estimation of the true
  functional gene content of subject genomes.
  - As a partial workaround, we recommend to check for paralogs of your query genes in your query genome before running BackBLAST so
    you are aware if there could be paralog-related issues. One easy way to check for paralogs is to run BLASTP for each of your query genes against
    the query genome and look for high-scoring hits. You can then watch for unexpected results when running BackBLAST. If you are not hitting
    real functional genes that you want to detect in the subject genomes due to having paralogs in the query genome, then you can mask out
    (i.e., delete) the paralogs of your query gene in the query genome file and re-run BackBLAST. Make sure you know what you are doing if
    you do this, though... in some use cases, you really want the paralogs to be there to filter out non-true hits to your query.
- If the __subject genomes__ contain paralogs (or multiple copies) of the query, then BackBLAST will not report this information. It will just
  show the best hit to your query gene and ignore the other hits. This is non-ideal if you want to know how many paralogs are in the subject genomes.
  - There is no built-in workaround for this in BackBLAST at the moment.

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
We hope that BackBLAST is helpful for you! If you use BackBLAST in your research, please cite the following paper,
which describes the initial version of BackBLAST:

> Bergstrand LH, Cardenas E, Holert J, Van Hamme JD, Mohn WW. Delineation of steroid-degrading microorganisms through
> comparative genomic analysis. __mBio__ 7:[10.1128/mbio.00166-16](https://doi.org/10.1128/mbio.00166-16).

# Appendix: full usage instructions
Help messages (e.g., from running `backblast -h`) are pasted below. These show some of the advanced options possible in
the command line (e.g., for `backblast auto`).

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
