BackBLAST-Gene-Cluster-Finder
==========================

A toolkit for finding and visualizing gene clusters across multiple bacterial genomes.

 - **BackBLAST.py** - Uses NCBI BLAST to search for gene clusters across multiple genomes. The script filtres out non-orthalagous genes from the search results by identifing and extracting only Bidirectional Best BLAST Hits using a graph-based algorithm. The algorithm is illistrated below:

![BackBLAST Algorithm] (https://raw.githubusercontent.com/LeeBergstrand/BackBLAST-Gene-Cluster-Finder/master/BackBLAST-Algorithm.gif)
