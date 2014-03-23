BackBLAST-Gene-Cluster-Finder
==========================

A toolkit for finding and visualizing gene clusters across multiple bacterial genomes.

 - **BackBLAST.py** - Uses NCBI BLAST to search for gene clusters across multiple genomes. The script filters out non-orthalagous genes from the search results by identifying and extracting only Bidirectional Best BLAST Hits using a graph-based algorithm. The algorithm is illustrated below:

![BackBLAST Algorithm] (https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/BackBLAST-Algorithm.gif)
