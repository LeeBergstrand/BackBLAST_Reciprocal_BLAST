BackBLAST_Reciprocal_BLAST
==========================
Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2019

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
This is a hacky, bleeding-edge install of BackBLAST, to be revised once we make a conda install.  
This method will only work in your single, active Bash session. One-time use.  
Keep checking in for updates to the code. Hoping to have a better solution soon.
```bash
# Download
# Before running this code, change directory into whatever folder you want to use the tool in
git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd BackBLAST_Reciprocal_BLAST
git checkout develop

# Install dependencies
conda env create -n backblast file="envs/conda_requirements.yaml"

# Add the repo scripts to your PATH temporarily in the current Bash session
PATH=${PATH}:${PWD}:${PWD}/scripts
```
Now you should be good to go! Run `BackBLAST -h` to get started.


# Usage
Rough notes on develop version for now.

## Recommended workflow
```bash
# Set up the run
BackBLAST.sh setup query.faa query_genome.faa subject_dir output_dir
# Then edit output_dir/config.yaml
# You can also edit output_dir/gene_metadata.tsv and output_dir/genome_metadata.tsv to make the plot look better

# Start the run
BackBLAST.sh run output_dir/config.yaml output_dir
# All done! You can iteratively refine the plot from here as you'd like.
```

## Speedy workflow
Gets the job done without any custom settings
```bash
BackBLAST.sh auto query.faa query_genome.faa subject_dir output_dir
```

## Test data
Try a test run from inside the repo with:
```bash
./BackBLAST.sh run ExampleData/Example_inputs/config.yaml ExampleData/Example_outputs --notemp
```

# Going deeper
For more thorough descriptions and information on usage please check the [**wiki**] (https://github.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/wiki) or look at the help information within `BackBLAST`. Enjoy!

