
library(igraph)

## load input data
setwd(R"(C:\Users\xiaod\Downloads\Rshiny\NGRN_Rshiny)") # set the working directory 
load("data_early_vs_late_demo.rdata") ## load "input_adj" and "a1"
save_graph_location <- r"(C:\Users\xiaod\Downloads\Rshiny\NGRN_Rshiny)"
getwd()

early_late_igraph <- function(adj_early,adj_late,prob_early,prob_late,ncRNA_num,ncRNA_node_color="pink",gene_node_color="slategray1",
                              gene_to_gene_color="tan1",ncRNA_to_gene_color="darkgray",
                              save_graph_location=getwd(),early_graph_name="early_edge_color_graph.jpeg",late_graph_name="late_edge_color_graph.jpeg"){
 
   # Take the element-wise maximum of adj_early and adj_late
  combined_adjacency <- pmax(adj_early, adj_late)
  # make sure prob_early and prob_late has weight numbers, + 0.2 for each matrix value
  prob_early <- prob_early+0.2
  prob_late <- prob_late+0.2
  
  ## test 1
  # Create igraph object from adj_early
  g_early <- graph_from_adjacency_matrix(adj_early, mode = "directed")
  
  # Create igraph object from adj_late
  g_late <- graph_from_adjacency_matrix(adj_late, mode = "directed")
  
  
  # Determine node type, color, and shape from combined_adjacency
  combined_node_type <- c(rep("ncRNA", ncRNA_num), rep("gene", ncol(adj_early) - ncRNA_num))
  combined_node_color <- ifelse(combined_node_type == "ncRNA", ncRNA_node_color, gene_node_color)
  combined_node_shape <- ifelse(combined_node_type == "ncRNA", "square", "circle")
  
  # Create a graph from the combined adjacency matrix
  g_combined <- graph_from_adjacency_matrix(combined_adjacency, mode = "directed")
  
  # Generate a layout based on the combined graph
  layout_combined <- layout_with_kk(g_combined)  # Kamada-Kawai layout for the combined graph
  
  
  # Determine edge colors based on source and target vertex types
  ecolor_early <- sapply(E(g_early), function(e) {
    edge_nodes <- ends(g_early, e, names = FALSE) ## this step specify the location in the location of a matrix
    source_vertex <- combined_node_type[edge_nodes[1]]
    target_vertex <- combined_node_type[edge_nodes[2]]
    if (source_vertex == "gene" && target_vertex == "gene") {
      return(gene_to_gene_color)
    } else if (source_vertex == "ncRNA" && target_vertex == "gene") {
      return(ncRNA_to_gene_color)
    } else {
      return("white")  # Default color for other types of connections
    }
  })
  
  ## see the first e in the code
  # edge_nodes <- ends(g_early, E(g_early)[1], names = FALSE)
  # source_vertex <- combined_node_type[edge_nodes[1]]
  # target_vertex <- combined_node_type[edge_nodes[2]]
  
  
  ecolor_late <- sapply(E(g_late), function(e) {
    edge_nodes <- ends(g_late, e, names = FALSE)
    source_vertex <- combined_node_type[edge_nodes[1]]
    target_vertex <- combined_node_type[edge_nodes[2]]
    if (source_vertex == "gene" && target_vertex == "gene") {
      return(gene_to_gene_color)
    } else if (source_vertex == "ncRNA" && target_vertex == "gene") {
      return(ncRNA_to_gene_color)
    } else {
      return("white")  # Default color for other types of connections
    }
  })
  
  # Adjust edge width based on weight matrices prob_early and prob_late
  ewidth_early <- sapply(E(g_early), function(e) {
    edge_nodes <- ends(g_early, e, names = FALSE)
    weight <- prob_early[edge_nodes[1], edge_nodes[2]]  # Get weight from prob_early matrix
    if (is.na(weight)) weight <- 1  # Default weight if NA (e.g., if no weight specified)
    return(weight)
  })
  
  ewidth_late <- sapply(E(g_late), function(e) {
    edge_nodes <- ends(g_late, e, names = FALSE)
    weight <- prob_late[edge_nodes[1], edge_nodes[2]]  # Get weight from prob_late matrix
    if (is.na(weight)) weight <- 1  # Default weight if NA (e.g., if no weight specified)
    return(weight)
  })
  
  setwd(save_graph_location)
  
  ## output is based on the size of the matrix
  n <- nrow(adj_early)
  edge.arrow.size <- ifelse(n <=45,1,
                            ifelse(n>45 & n<=75, 0.5,
                                   ifelse(n>75, 0.4,0.2)))
  edge.width_num <- ifelse(n <=45,4,
                           ifelse(n>45 , 2, 1))
  vertex.size <- ifelse(n <=45,12,
                     ifelse(n>45 & n<=75, 7,
                            ifelse(n>75, 5,5)))
  vertex.label.cex <- ifelse(n <=45,1,
                             ifelse(n>45 & n<=75, 0.7,
                                    ifelse(n>75, 0.7,0.7)))
  
  
  # Plot graph from adj_early with edge colors based on combined_adjacency
  # tiff(early_graph_name, units = "in", width = 13, height = 10, res = 100) 
  plot(g_early, vertex.color = combined_node_color, vertex.shape = combined_node_shape,
       vertex.label.color = "black", layout = layout_combined, edge.arrow.size = edge.arrow.size,
       edge.color = ecolor_early, edge.width = edge.width_num*ewidth_early, vertex.size = vertex.size,vertex.label.cex = vertex.label.cex)
  dev.off()
  
  # Plot graph from adj_late with edge colors based on combined_adjacency
  # tiff(late_graph_name, units = "in", width = 13, height = 10, res = 100) 
  plot(g_late, vertex.color = combined_node_color, vertex.shape = combined_node_shape, 
       vertex.label.color = "black", layout = layout_combined, edge.arrow.size = edge.arrow.size, 
       edge.color = ecolor_late, edge.width = edge.width_num*ewidth_late, vertex.size = vertex.size,vertex.label.cex = vertex.label.cex) 
  # dev.off()
}

## test the function
p=6 # determine how many ncRNA in the matrix first # first 6 rows
early_late_igraph(adj_early,adj_late,prob_early,prob_late,ncRNA_num=p,early_graph_name="early_edge_color_graph2.jpeg",late_graph_name="late_edge_color_graph2.jpeg")
