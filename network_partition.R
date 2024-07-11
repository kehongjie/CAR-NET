## NOTE: This file contains the codes to partition a large network into 
##        some smaller sub-network (also refereed as modules). 

## required package:
library(igraph)

## INPUT:
## Suppose you have a adjacency matrix called "adj" of dimension (p1+q1)*(p1+q1),
##    where p1 is the number of ncRNAs in the network and q1 is the number of genes
##    in the network. 
## And you also have a posterior matrix called "mat_post" of dimension (p1+q1)*(p1+q1)

partition <- function(adj, mat_post, p1) {
  ## network partition
  g1 <- graph_from_adjacency_matrix(adj, mode="undirected")
  set.seed(1234)
  fit_lou <- cluster_louvain(g1, resolution=1)
  # sort(table(fit_lou$membership), decreasing = T)
  
  ## Now extract the network partition results
  ## fit_lou$membership gives the module assignment infor; for example, if its first 
  ##    element is 5, that means first node is assigned to fifth modules; 
  ##    length(fit_lou$membership) will equal to nrow(adj)
  ## grp_idx tells you the names of those largest modules, in descending order; 
  ##    for example, if its first element is 5, that means the fifth module is the 
  ##    largest module (has the most nodes).
  grp_idx <- as.numeric(names(sort(table(fit_lou$membership), 
                                   decreasing = T))) ## group name in descending order
  grp_idx <- grp_idx[as.numeric(sort(table(fit_lou$membership), 
                                     decreasing = T))>=10] ## only keep modules with size larger than 10
  
  return(list(adj=adj, mat_post=mat_post, fit=fit_lou$membership, p=p1, grp_idx=grp_idx))
}
