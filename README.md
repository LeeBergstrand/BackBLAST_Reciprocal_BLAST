BackBLAST-Gene-Cluster-Finder
==========================

A toolkit for finding and visualizing gene clusters across multiple bacterial genomes.

- **BackBLAST.py** - A script that uses NCBI BLAST to search for gene clusters across multiple genomes. Non-orthalagous genes are filtred out by identifying and extracting only Bidirectional Best BLAST Hits using a graph-based algorithm. The algorithm is illustrated below:
![BackBLAST Algorithm] (https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/BackBLAST-Algorithm.gif)
- **Visualization** - This repository also includes tools for visualizing the results of **BackBLAST.py** in the form of a R heatmap. These tools are located in the visualization [sub-directory](https://github.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/tree/master/Visualization).

For more thorough descriptions and information on usage please check the [**wiki!**] (https://github.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/wiki)
