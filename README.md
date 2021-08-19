BackBLAST_Reciprocal_BLAST
==========================
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3465955.svg)](https://doi.org/10.5281/zenodo.3465955)

Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2021

This repository contains a reciprocal BLAST program for filtering down BLAST results to best bidirectional hits. It also contains a toolkit for finding and visualizing BLAST hits for gene clusters within multiple bacterial genomes.

# Repo contents
(To be updated once BackBLAST2 is complete)

- **BackBLAST_core.py** - A script that uses NCBI BLAST to search for gene clusters within a within a bacterial genome genome. Non-orthalagous genes are filtred out by identifying and extracting only Bidirectional BLAST Hits using a graph-based algorithm. The algorithm is illustrated below:

![BackBLAST Algorithm](https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/Media/BackBLAST-Algorithm.gif)

- **Visualization** - This repository also includes tools for visualizing the results from **BackBLAST_core.py** in the form of a R heatmap. Here is an example:

![Example Results](https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/Media/ExampleResults.jpeg)


# Installation
Temporary instructions while BackBLAST2 is still under development, to install the development version  

## Dependencies
- Linux operating system (e.g., Ubuntu)
- miniconda (2 or 3) or Anaconda
- Workflow is pretty light on RAM, CPU, and storage space, so most machines should be able to handle BackBLAST without issue. The only exception is if you create genome trees within the pipeline, in which case you'll need a fair amount of CPU and time to calculate large trees.

## Instructions
Takes a few steps -- to be revised once we make a conda install.
```bash
# Download the repo
cd /tmp
git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd BackBLAST_Reciprocal_BLAST
git checkout develop # optionally go to a specific branch

# Create the conda env based on the YAML file in the repo
# It is recommended that you run this command using mamba instead of conda - conda might fail during install.
mamba env create -n backblast --file=envs/conda_requirements.yaml

# Copy the key repo contents into a conda share folder
conda activate backblast
mkdir -p ${CONDA_PREFIX}/share/BackBLAST
cp -r * ${CONDA_PREFIX}/share/BackBLAST

# Remove the original repo
cd ..
rm -rf BackBLAST_Reciprocal_BLAST

# Add instructions to export the BackBLAST folder to the PATH when the repo activates
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

if [[ ! -f ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh ]]; then
  echo '#!/bin/sh' > ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
fi
echo "export PATH=\${PATH}:${CONDA_PREFIX}/share/BackBLAST:${CONDA_PREFIX}/share/BackBLAST/scripts" \
  >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
chmod 755 ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# Re-activate the repo to apply the changes
conda deactivate
conda activate backblast
```
Now you should be good to go! Run `BackBLAST -h` to get started.


# Usage
Rough notes on develop version for now.

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
backblast auto query.faa query_genome.faa subject_dir output_dir
```

## Test data
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

# Going deeper
For more thorough descriptions and information on usage please check the [**wiki**] (https://github.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/wiki) or look at the help information within `BackBLAST`. Enjoy!

