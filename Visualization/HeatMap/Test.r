library(ape)

tree = read.tree("/Users/lee/Dropbox/RandD/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/SteriodBacteria16SAligned2.nwk")

print("Please select a root.")
plot(tree)

tree = root(tree, resolve.root = TRUE, interactive = TRUE)

print("Thanks, continuing...")
plot(tree)

tree$edge.length[which(tree$edge.length == 0)] = 0.00001
tree <- chronopl(tree, lambda = 0.1, tol = 0)
plot(tree)

dendrogram = as.dendrogram(as.hclust.phylo(tree))
plot(dendrogram)

cladeOrder = order.dendrogram(dendrogram)
cladeName =  labels(dendrogram)
cladePosition = data.frame(cladeName, cladeOrder)
print(cladePosition)
cladePosition = cladePosition[order(cladePosition$cladeOrder),]
print(cladePosition)