BackBLAST_Reciprocal_BLAST
==========================
Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019

This repository contains a reciprocal BLAST program for filtering down BLAST results to best bidirectional hits. It also contains a toolkit for finding and visualizing BLAST hits for gene clusters within multiple bacterial genomes.

# Repo contents
(To be updated once BackBLAST2 is complete)

- **BackBLAST.py** - A script that uses NCBI BLAST to search for gene clusters within a within a bacterial genome genome. Non-orthalagous genes are filtred out by identifying and extracting only Bidirectional BLAST Hits using a graph-based algorithm. The algorithm is illustrated below:

![BackBLAST Algorithm](https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/Media/BackBLAST-Algorithm.gif)

- **Visualization** - This repository also includes tools for visualizing the results from **BackBLAST.py** in the form of a R heatmap. Here is an example:

![Example Results](https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/Media/ExampleResults.jpeg)


# Installation
Temporary instructions while BackBLAST2 is still under development, to install the development version  

## Dependencies
- Linux operating system (e.g., Ubuntu)
- miniconda (2 or 3) or Anaconda
- Workflow is pretty light on RAM, CPU, and storage space, so most machines should be able to handle BackBLAST without issue. The only exception is if you create genome trees within the pipeline, in which case you'll need a fair amount of CPU and time to calculate large trees.

## Instructions
```bash
# Download
cd /tmp
git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd /tmp/BackBLAST_Reciprocal_BLAST
git checkout develop

# Install dependencies
conda create -n backblast -c bioconda -c conda-forge snakemake=5.5.4

# Install main scripts
conda activate backblast
cp -r BackBLAST.sh BackBLAST.py Graph.py Snakefile template_config.yaml envs ${CONDA_PREFIX}/bin
cd Visualization
cp RemoveDuplicates.sh CreateBlankResults.py CombineBlastTables.R generate_BackBLAST_heatmap.R ${CONDA_PREFIX}/bin

# Clean up
cd /tmp && rm -r BackBLAST_Reciprocal_BLAST
```
Now you should be good to go! Run `BackBLAST.sh -h` to get started.


# Usage
Rough notes on develop version for now

## Recommended workflow
```bash
# Set up the run
BackBLAST.sh setup query.faa query_genome.faa subject_dir/ output_dir/
# Then edit output_dir/config.yaml
# You can also edit output_dir/gene_metadata.tsv and output_dir/genome_metadata.tsv to make the plot look better

# Start the run
BackBLAST.sh run output_dir/
# All done! You can iteratively refine the plot from here as you'd like.
```

## Speedy workflow
Gets the job done without any custom settings
```bash
BackBLAST.sh auto query.faa query_genome.faa subject_dir/ output_dir/
```

# Going deeper
For more thorough descriptions and information on usage please check the [**wiki**] (https://github.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/wiki) or look at the help information within `BackBLAST.sh`. Enjoy!

