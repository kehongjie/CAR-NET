
######################################################################
## adj generator for single matrix (adj and prop)
######################################################################

## FUNCTION: to simulate a (ncRNA-gene) adjacency matrix, given the number of 
####  ncRNAs, number of genes, and number of edges (both types). 

## p: number of ncRNAs
## q: number of genes
## l1: number of level-1 edges (output could be larger than this number)
## l2: number of level-2 edges (output could be larger than this number)

adj.generator <- function(p, q, l1, l2) {
  out <- matrix(0, p+q, p+q) 
  colnames(out) <- rownames(out) <- c(paste("ncRNA", 1:p, sep=""),
                                      paste("gene", 1:q, sep=""))
  
  if (p>q) {stop("There should be less ncRNAs than genes!!")}
  if (l1<=p) {stop("There should be more level-1 edges!!")}
  if (l1>(p*q)) {stop("Too many level-1 edges!!")}
  if (l2>(q*(q-1))) {stop("Too many level-2 edges!!")}
  
  tmp <- matrix(0, p, q)
  diag(tmp) <- 1 # make sure all the ncRNAs have at least one edge
  pos <- which((upper.tri(tmp) | lower.tri(tmp)) ==T)
  tmp[sample(pos, l1-p)] <- 1
  out[1:p, (p+1):(p+q)] <- tmp # level-1 edges
  
  out_sub <- out[(p+1):(p+q), (p+1):(p+q)]
  pos <- which((upper.tri(out_sub) | lower.tri(out_sub)) ==T)
  out[(p+1):(p+q), (p+1):(p+q)][sample(pos, l2)] <- 1 # level-2 edges 
  
  pos_check <- union(setdiff(which(colSums(out)==0), 1:p), which(rowSums(out)==0))
  if (length(pos_check)>0) {
    for (idx in pos_check) {
      out[sample(setdiff(1:(p+q),idx), 1), idx] <- 1 
    }
  } # make sure all the genes have at least one edge
  
  return(out)
}

######################################################################
## function for single matrix
######################################################################


single_igraph <- function(adj_single,prob_single,ncRNA_num,ncRNA_node_color="#ffc000",gene_node_color="#01b0f0",
                                        gene_to_gene_color="darkgray",ncRNA_to_gene_color="darkgray",
                                        save_graph_location=getwd(),single_graph_name="single_edge_color_graph.jpeg"){
  
  # make sure prob_single has weight numbers, + 0.2 for each matrix value
  prob_single <- prob_single+0.2
  
  
  library(igraph)
  # Create igraph object from adj_single
  g_single <- graph_from_adjacency_matrix(adj_single, mode = "directed")
  
  # Determine node type, color, and shape from adj_single
  node_type <- c(rep("ncRNA", ncRNA_num), rep("gene", ncol(adj_single) - ncRNA_num))
  node_color <- ifelse(node_type == "ncRNA", ncRNA_node_color, gene_node_color)
  # node_shape <- ifelse(node_type == "ncRNA", "square", "circle")
  node_shape <- rep("circle", length(node_type)) ## all circle shape 
  
  # Determine edge colors based on source and target vertex types
  ecolor_single <- sapply(E(g_single), function(e) {
    edge_nodes <- ends(g_single, e, names = FALSE) ## this step specify the location in the location of a matrix
    source_vertex <- node_type[edge_nodes[1]]
    target_vertex <- node_type[edge_nodes[2]]
    if (source_vertex == "gene" && target_vertex == "gene") {
      return(gene_to_gene_color)
    } else if (source_vertex == "ncRNA" && target_vertex == "gene") {
      return(ncRNA_to_gene_color)
    } else {
      return("white")  # Default color for other types of connections
    }
  })
  
  
  # Adjust edge width based on weight matrices prob_single
  ewidth_single <- sapply(E(g_single), function(e) {
    edge_nodes <- ends(g_single, e, names = FALSE)
    weight <- prob_single[edge_nodes[1], edge_nodes[2]]  # Get weight from prob_early matrix
    if (is.na(weight)) weight <- 1  # Default weight if NA (e.g., if no weight specified)
    return(weight)
  })
  
  ## output is based on the size of the matrix
  n <- nrow(adj_single)
  edge.arrow.size <- ifelse(n <=10,1,
                            ifelse(n>10 & n<=40, 0.8,
                                   ifelse(n>40 & n<=75, 0.5, 0.4)))
  edge.width_num <- ifelse(n <=45,4,
                           ifelse(n>45 , 2, 1))
  vertex.size <- ifelse(n <=45,12,
                        ifelse(n>45 & n<=75, 7,
                               ifelse(n>75, 5,5)))
  vertex.label.cex <- ifelse(n <=45,1,
                             ifelse(n>45 & n<=75, 0.7,
                                    ifelse(n>75, 0.7,0.7)))
  
  setwd(save_graph_location)
  # Plot graph from adj_early with edge colors based on combined_adjacency
  # tiff(single_graph_name, units = "in", width = 10, height = 10, res = 150)
  plot(g_single, vertex.color = node_color, vertex.shape = node_shape, layout = layout_with_fr(g_single),
       vertex.label.color = "black", layout = layout_with_kk, edge.arrow.size = edge.arrow.size, 
       edge.color = ecolor_single, edge.width = edge.width_num*ewidth_single, 
       vertex.size = vertex.size,vertex.label.cex = vertex.label.cex) 
  
  # f <- factor(levels = c("ncRNA Node", "gene Node", "gene to gene", "ncRNA to gene"), ordered = TRUE)
  # vcols <- c("pink", "slategray1", "tan1", "darkgray")
  f <- factor(levels = c("ncRNA node", "gene node"), ordered = TRUE)
  vcols <- c(ncRNA_node_color, gene_node_color)
  legend("topleft", legend = levels(f), pch = 19, col = vcols, bty = "n")
  
  # dev.off()
}

######################################################################
## test: output graphs
######################################################################


## example: 0-45
# p <- 10
# q <- 25
# adj <- adj.generator(p=p, q=q, l1=12, l2=15)
# sum(adj) # check number of edges 
# sum(adj[1:p, (p+1):(p+q)]) # check number of level-1 edges
# sum(adj[(p+1):(p+q), (p+1):(p+q)]) # check number of level-2 edges
# 
# single_igraph(adj_single=adj, prob_single=matrix(0.5, p+q, p+q), 
#                          ncRNA_num=p,
#                          single_graph_name="single_edge_color_graph_0613_0-45.jpeg") # call your visulaization function
# 
# ## example: 46-75
# p <- 30
# q <- 40
# adj <- adj.generator(p=p, q=q, l1=35, l2=20)
# sum(adj) # check number of edges 
# sum(adj[1:p, (p+1):(p+q)]) # check number of level-1 edges
# sum(adj[(p+1):(p+q), (p+1):(p+q)]) # check number of level-2 edges
# 
# single_igraph(adj_single=adj, prob_single=matrix(0.5, p+q, p+q), 
#                          ncRNA_num=p,
#                          single_graph_name="single_edge_color_graph_0613_46-75.jpeg") # call your visulaization function
# 
# 
# ## example: 76 and above
# p <- 40
# q <- 40
# adj <- adj.generator(p=p, q=q, l1=45, l2=45)
# sum(adj) # check number of edges 
# sum(adj[1:p, (p+1):(p+q)]) # check number of level-1 edges
# sum(adj[(p+1):(p+q), (p+1):(p+q)]) # check number of level-2 edges
# 
# single_igraph(adj_single=adj, prob_single=matrix(0.5, p+q, p+q), 
#                             ncRNA_num=p,
#                             single_graph_name="single_edge_color_graph_0613_76-above.jpeg") # call your visulaization function
