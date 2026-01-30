all the genes that are present in a tissue

breast tissue today
3 billion genes total
2000 genes in the readout

if someone had cancer vs no cancer. what is the 

GTEx_network: contains correlated genes and gene groups, which tissue they express in, and the weight of their correlation.
GTEx_module: each gene symbol found in GTEx_network is catalogued here along with its tissue, module, layer.

its still unclear to me what REACTOME is and why the geneIDs are so long but for now we will focus on GTEx_network and GTEx_module.

Paper:
gene co-expression network analysis has been shown effective in identifying functional co-expressed gene modules associated with complex human diseases. 

Existing techniques to construct co-expression networks require some critical prior information such as predefined number of clusters, numerical thresholds for defining co-expression/interactions
or they dont naturally reproduce the hallmarks of complex systems such as "scale-free degree distribution fo small worldness".[1]

i. e. they dont naturally reproduce the hallmarks of any complex systems. like that in complex systems, nodes arent directly connected to every other node, but is a short number of hops from most node.

in the past a graph filtering technique called planar maximally filtered graph(PMFG) has been applied to real world datasets, but it does not fit in this use case. due to its complexity O(|V|^3)

this paper provides a new framework called multiscale embedded gene co-expression network analysis.(MEGENA) it works by 
1. introducing quality control of co-expression similarities
2. parallelizing embedded network construction
3. developing a novel clustering technique to identify multi-scale clustering structures in Planar Filtered Networks. 




[1]small-worldness is a type of mathematical graph in which most nodes are not neighbors of one another, but the neighbors of any given node are likely to be neighbors of each other and most nodes can be reached from every other node via a small number of hops.

df.apply sounds like a weird but helpful function

ok so we have 50 tissues. Look at graph-metrics.csv for more info and look at tissue_networks/ for the networks themselves.


TODO: 
- look at MEGENA  Paper
- find out more about the bounds of visualizing graphs