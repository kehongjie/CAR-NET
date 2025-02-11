saved_data_server <- function(input, output, session, ImProxy) {
  
  ns <- NS("saved_data")
  
  ##########################
  # Reactive Values        #
  ##########################
  DB <- reactiveValues(
    MergedDB=MergedDB.load(db), 
    MergedSpecies=MergedSpecies.load(db), 
    MergedStudyNames=MergedStudyNames.load(db), 
    pathway.list=NULL,
    ACS_ADS_global=NULL,
    ACS_ADS_pathway=NULL
  )
  
  adj_mat <- reactiveVal(value = NULL)  # final adjacency matrix
  post_prob <- reactiveVal()  # posterior probability matrix
  p <- reactiveVal()  # number of ncRNA
  q <- reactiveVal()  # number of genes
  
  shared_var <- reactiveValues()
  
  ##########################
  # Observers              #
  ##########################
  
  # Bayesian Network Generation
  source("./internal_functions/bn_main.R")
  source("./internal_functions/gui_filter.R")
  source("./internal_functions/single-adj_different_size.r")
  source("./internal_functions/network_partition.R")
  source("./internal_functions/run_pathway.R")
  
  ########## run BN and overview ##########
  observeEvent(input$run_algo, {
    wait(session, "Generating Bayesian Network, may take a while")
    path_old <- getwd()
    try({
      filePath <- ImProxy$file1
      fileText <- read.csv(filePath$datapath, check.names = FALSE)
      df <- t(fileText)
      colnames(df) <- df[1,]
      df <- df[-1,]
      df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
      X <- as.matrix(internal_filter(df, ImProxy$cutoff_ncRNA)[1][[1]])
      #
      filePath2 <- ImProxy$file2
      fileText2 <- read.csv(filePath2$datapath, check.names = FALSE)
      df <- t(fileText2)
      colnames(df) <- df[1,]
      df <- df[-1,]
      df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
      Y <- as.matrix(internal_filter(df, ImProxy$cutoff_gene)[1][[1]])
      
      shared_var$X <- X
      shared_var$Y <- Y
      
      # x <- read.csv("./data/TCGA_KIRP_early_lncRNA_top5per.csv", header=T)
      # rownames(x) <- x[,1]
      # x <- x[,-1]
      # X <- t(x)
      # X <- as.matrix(internal_filter(X, 3)[1][[1]])
      # y <- read.csv("./data/TCGA_KIRP_early_gene_top5per.csv", header=T)
      # rownames(y) <- y[,1]
      # y <- y[,-1]
      # Y <- t(y)
      # Y <- as.matrix(internal_filter(Y, 11)[1][[1]])

      fit <- bn.main(X=X, Y=Y, alpha=1e-5)
      # load("./bn.RData")
      
      # output$text <- renderText({ paste("After Stage 1:\n",
      #   "Left with ", fit$s1_lvl1, " level-1 edges\n",
      #   "Left with ", fit$s1_lvl2, " level-2 edges\n",
      #   "Involving ", fit$s1_rna, " ncRNAs and ", fit$s1_gene, " genes\n",
      #   "After Stage 2:\n",
      #   "Left with ", fit$s2_lvl1, " level-1 edges\n", 
      #   "Left with ", fit$s2_lvl2, " level-2 edges\n",
      #   "Involving ", fit$s2_rna, " ncRNAs and ", fit$s2_gene, " genes", sep="") })
      
      ## output for panel 1
      output$text <- renderText({paste("Final network inlcudes ", 
                                       fit$s2_rna, " ncRNAs and ", fit$s2_gene, " genes;\n",
                                       fit$s2_lvl1, " ncRNA-->ncRNA edges and ",
                                       fit$s2_lvl2, " gene-->gene edges.",
                                       sep="") })
      
      output$table_pair <- renderTable(fit$mat_pair, striped=TRUE)
      
      
      # rownames(fit$adj) <- c(ImProxy$names_rna, ImProxy$names_gene)
      # colnames(fit$adj) <- c(ImProxy$names_rna, ImProxy$names_gene)
      
      pos_innet <- (rowSums(fit$adj)>0 | colSums(fit$adj)>0) ## position of in-network nodes
      small <- fit$adj[pos_innet, pos_innet]
      dag_sub <- fit$adj[pos_innet, pos_innet] ## remove rows/columns that are all zero
      adj_mat(dag_sub) ## dimension p2*q2
      post_prob(fit$post.prob[pos_innet, pos_innet])
      p1 <- fit$s1_rna
      q1 <- fit$s1_gene
      p2 <- sum(pos_innet[1:p1])
      q2 <- sum(pos_innet[(p1+1):(p1+q1)])
      
      shared_var$pos_innet <- pos_innet
      shared_var$p2 <- p2
      shared_var$q2 <- q2
      
      save(fit, file="bn.RData")
      print("Bayesian Network saved as 'bn.RData'.")
      
  
      
      sendSuccessMessage(session, "Bayesian Network generated.")
    }, session)
    
    observe({
      updateTabsetPanel(session, "tabSelect",
                        selected = "panel1")
    })
    done(session)
  })
  
  ########## pathway analysis for all genes in the network ##########
  observeEvent(input$pathway_all, {
    p2 <- shared_var$p2
    q2 <- shared_var$q2
    
    gene_list <- colnames(adj_mat())[(p2+1):(p2+q2)]
    bg <- colnames(shared_var$Y)
    print(paste(length(gene_list), "genes in the network"))
    print(paste(length(bg), "genes in the background"))
    
    output$path_plot_all <- renderPlot({
      run.path(gene.list=gene_list, background=bg)
    })
    
    observe({
      updateTabsetPanel(session, "tabSelect",
                        selected = "panel1")
    })
    
  }) 
  
  
  ########## modulization ##########
  observeEvent(input$butt_mod, {
    
    mat_post_sub <- post_prob()
    p2 <- shared_var$p2
    q2 <- shared_var$q2
    print(paste("p2=", p2, sep=""))
    print(paste("q2=", q2, sep=""))
    
    # result <- partition(adj_mat(), post_prob()[pos_innet, pos_innet], p())
    ## network partition
    g1 <- graph_from_adjacency_matrix(adj_mat(), mode="undirected")
    set.seed(1234)
    fit_lou <- cluster_louvain(g1, resolution=1)
    grp_idx <- as.numeric(names(sort(table(fit_lou$membership), 
                                     decreasing = T))) ## group name in descending order
    grp_idx <- grp_idx[as.numeric(sort(table(fit_lou$membership), 
                                       decreasing = T))>=5] ## only keep modules with size larger than 10
    shared_var$fit_lou <- fit_lou
    shared_var$grp_idx <- grp_idx
    print(paste("length of membership is", length(fit_lou$membership)))
    print(fit_lou$membership)
    print(grp_idx)
    
    ## only visualize a subset of modules 
    m <- tmp <- 0
    mod_size <- as.numeric(sort(table(fit_lou$membership), decreasing = T))
    while (tmp<100 & m<length(mod_size)) {
      m <- m+1
      tmp <- sum(mod_size[1:m])
    } ## take up to 100 nodes
    if(tmp>100) {m <- m-1}
    print(paste("m=",m,sep="")) 
    print(mod_size)
    
    ## reorder nodes within modules (later)
    # ncord <- geord <- ncsep <- gesep <- NULL
    # for (i in 1:m) {
    #   pos_nc <- intersect(which(fit_lou$membership==grp_idx[i]), 1:p2)
    #   pos_ge <- intersect(which(fit_lou$membership==grp_idx[i]), (p2+1):(p2+q2))
    #   
    #   post_mod <- mat_post_sub[pos_nc, pos_ge]
    # 
    #   fitp <- pheatmap(post_mod, cluster_rows = T, cluster_cols = T, fontsize_row=5, fontsize_col=5,
    #                    color = colorRampPalette(c("white","red4"))(50), silent=T)
    #   
    #   ncord <- c(ncord, rownames(post_mod)[fitp$tree_row$order])
    #   geord <- c(geord, colnames(post_mod)[fitp$tree_col$order])
    #   ncsep <- c(ncsep, length(ncord))
    #   gesep <- c(gesep, length(geord))
    # }
    
    ## not reordering 
    ncord <- geord <- ncsep <- gesep <- NULL
    mat_node <- NULL
    for (i in 1:m) {
      pos_nc <- intersect(which(fit_lou$membership==grp_idx[i]), 1:p2)
      pos_ge <- intersect(which(fit_lou$membership==grp_idx[i]), (p2+1):(p2+q2))
      
      ncord <- c(ncord, rownames(mat_post_sub)[pos_nc])
      geord <- c(geord, colnames(mat_post_sub)[pos_ge])
      ncsep <- c(ncsep, length(ncord))
      gesep <- c(gesep, length(geord))
      
      ## for node names table
      tmp1 <- rownames(mat_post_sub)[pos_nc]
      tmp2 <- colnames(mat_post_sub)[pos_ge]
      new_mat <- cbind(rep(paste("Module",i), length(tmp1)+length(tmp2)), 
                       c(tmp1, tmp2),
                       c(rep("ncRNA",length(tmp1)), rep("gene",length(tmp2))))
      mat_node <- rbind(mat_node, new_mat)
      if (i<m) {mat_node <- rbind(mat_node, c(rep(" ",3)), 
                                  c(rep(" ",3)))} ## for separating modules
    }
    colnames(mat_node) <- c("Module", "Node name", "Node type")

    
    post_mod <- mat_post_sub[match(ncord, colnames(mat_post_sub)), 
                             match(geord, colnames(mat_post_sub))]
    
    output$heatmap <- renderPlot({
      pheatmap(post_mod, cluster_rows = F, cluster_cols = F, fontsize_row=8, fontsize_col=8,
               color = colorRampPalette(c("#f4f5f4", "#b2182b"))(50), border_color="black",
               gaps_row=ncsep, gaps_col=gesep, angle_col=315) # fontsize 3-6
      
    })
    
    output$table_node <- renderTable(mat_node, striped=TRUE)
    
    observe({
      updateTabsetPanel(session, "tabSelect",
                        selected = "panel2")
    })
    
    
  })
  
  
  
  ########## module visulization (igraph) ##########
  observeEvent(input$butt_mod, {
    # wait(session, "Generating Modules")
    
    ## update drop-down input
    # try({
      # pos_innet <- (rowSums(adj_mat())>0 | colSums(adj_mat())>0)
      # result <- partition(adj_mat(), post_prob()[pos_innet, pos_innet], p())
    fit_lou <- shared_var$fit_lou 
    grp_idx <- shared_var$grp_idx
    mat_post_sub <- post_prob()
    dag_sub <- adj_mat()
    p2 <- shared_var$p2
    q2 <- shared_var$q2
    
    l <- as.list(seq(1,length(grp_idx)))
    observe({
      updateSelectInput(session, "sel_mod",
                        label = "Select a module to display",
                        choices = l,
                        selected = l[1])
      })
    
    output$network_visual <- renderPlot({
      # # A temp file to save the output.
      # # This file will be removed later by renderImage
      # outfile <- tempfile(fileext = '.png')
      # # Generate the PNG
      # png(outfile, width = 800, height = 600)
      
      idx <- grp_idx[as.numeric(input$sel_mod)]
      pos_node <- which(fit_lou$membership==idx) ## position of nodes in this module
      adj_single <- dag_sub[pos_node, pos_node]
      prob_single <- mat_post_sub[pos_node, pos_node]
      a1 <- length(intersect(pos_node, 1:p2)) ## number of ncRNAs in this module
      
      ## NOTE: modify the visualization function accordingly
      single_igraph(adj_single=adj_single, prob_single=prob_single, ncRNA_num=a1,
                    single_graph_name="modify_the_graph_name.jpeg")
      # dev.off()
      
      # Return a list containing the filename
      # list(src = outfile,
      #      contentType = 'image/png',
      #      width = 800,
      #      height = 600,
      #      alt = "This is alternate text")
    })#, deleteFile = TRUE)
    
    
    observe({
      updateTabsetPanel(session, "tabSelect",
                        selected = "panel3")
    })
    done(session)
  })
  
  # end
  
}
