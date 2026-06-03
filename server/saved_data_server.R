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
      if (is.null(ImProxy$file1) || is.null(ImProxy$file2) ||
          is.null(ImProxy$filtered_ncRNA_data) || is.null(ImProxy$filtered_gene_data)) {
        stop("Select TCGA_KIRP example or upload both ncRNA and gene expression CSVs before running CAR-NET.")
      }
      filePath <- ImProxy$file1
      fileText <- read.csv(filePath$datapath, check.names = FALSE)
      df <- t(fileText)
      colnames(df) <- df[1,]
      df <- df[-1,]
      df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
      #browser()
      X <- df
      #X <- t(ImProxy$filtered_ncRNA_data)
      
      #
      filePath2 <- ImProxy$file2
      fileText2 <- read.csv(filePath2$datapath, check.names = FALSE)
      df <- t(fileText2)
      colnames(df) <- df[1,]
      df <- df[-1,]
      df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
      Y <- t(ImProxy$filtered_gene_data)
      #Y <- df 
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
  
  ########## modulization ##########
  observeEvent(input$butt_mod, {
    mat_post_sub <- post_prob()
    dag_sub <- adj_mat()
    p2 <- shared_var$p2
    q2 <- shared_var$q2
    min_module_size <- as.numeric(input$mod_size)
    if (!is.finite(min_module_size) || min_module_size < 1) {
      min_module_size <- 1
    }
    min_module_size <- ceiling(min_module_size)
    
    if (is.null(mat_post_sub) || is.null(dag_sub) ||
        is.null(p2) || is.null(q2) || nrow(dag_sub) == 0) {
      output$heatmap <- renderPlot({
        plot.new()
        text(0.5, 0.5, "Run CAR-NET before detecting modules.")
      })
      output$table_node <- renderTable(
        data.frame(Message = "Run CAR-NET before detecting modules."),
        striped = TRUE
      )
      output$network_visual <- renderPlot({
        plot.new()
        text(0.5, 0.5, "Run CAR-NET before detecting modules.")
      })
      updateSelectInput(session, "sel_mod", choices = c("No modules available" = ""), selected = "")
      updateTabsetPanel(session, "tabSelect", selected = "panel2")
      return(NULL)
    }
    
    nc_positions <- if (p2 > 0) seq_len(p2) else integer(0)
    gene_positions <- if (q2 > 0) (p2 + 1):(p2 + q2) else integer(0)
    
    # result <- partition(adj_mat(), post_prob()[pos_innet, pos_innet], p())
    ## network partition
    g1 <- graph_from_adjacency_matrix(dag_sub, mode="max")
    set.seed(1234)
    fit_lou <- cluster_louvain(g1, resolution=1)
    module_sizes <- sort(table(fit_lou$membership), decreasing = TRUE)
    grp_idx <- as.numeric(names(module_sizes)[as.numeric(module_sizes) >= min_module_size])
    has_ncRNA_and_gene <- vapply(grp_idx, function(idx) {
      pos_node <- which(fit_lou$membership == idx)
      length(intersect(pos_node, nc_positions)) > 0 &&
        length(intersect(pos_node, gene_positions)) > 0
    }, logical(1))
    grp_idx <- grp_idx[has_ncRNA_and_gene]
    shared_var$fit_lou <- fit_lou
    shared_var$grp_idx <- grp_idx
    
    if (length(grp_idx) == 0) {
      msg <- paste0(
        "No modules with at least ", min_module_size,
        " nodes and both ncRNA/gene nodes were found. Try a smaller minimum module size."
      )
      output$heatmap <- renderPlot({
        plot.new()
        text(0.5, 0.5, msg)
      })
      output$table_node <- renderTable(data.frame(Message = msg), striped = TRUE)
      output$network_visual <- renderPlot({
        plot.new()
        text(0.5, 0.5, msg)
      })
      updateSelectInput(session, "sel_mod", choices = c("No modules available" = ""), selected = "")
      updateTabsetPanel(session, "tabSelect", selected = "panel2")
      return(NULL)
    }
    
    ## only visualize a subset of modules 
    kept_module_sizes <- as.numeric(module_sizes[match(grp_idx, names(module_sizes))])
    cumulative_size <- cumsum(kept_module_sizes)
    m <- max(which(cumulative_size <= 100), na.rm = TRUE)
    if (!is.finite(m) || m < 1) {
      m <- 1
    }
    m <- min(m, length(grp_idx))
    display_grp_idx <- grp_idx[seq_len(m)]
    
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
    for (i in seq_along(display_grp_idx)) {
      pos_nc <- intersect(which(fit_lou$membership == display_grp_idx[i]), nc_positions)
      pos_ge <- intersect(which(fit_lou$membership == display_grp_idx[i]), gene_positions)
      
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

    post_mod <- mat_post_sub[match(ncord, rownames(mat_post_sub)), 
                             match(geord, colnames(mat_post_sub))]
    
    output$heatmap <- renderPlot({
      validate(need(nrow(post_mod) > 0 && ncol(post_mod) > 0,
                    "No module heatmap is available for the selected minimum module size."))
      pheatmap(post_mod, cluster_rows = F, cluster_cols = F, fontsize_row=8, fontsize_col=8,
               color = colorRampPalette(c("#f4f5f4", "#b2182b"))(50), border_color="black",
               gaps_row=head(ncsep, -1), gaps_col=head(gesep, -1), angle_col=315) # fontsize 3-6
      
    })
    
    output$table_node <- renderTable(mat_node, striped=TRUE)
    
    updateTabsetPanel(session, "tabSelect", selected = "panel2")
    
    
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
    
    if (is.null(fit_lou) || is.null(grp_idx) || length(grp_idx) == 0) {
      updateSelectInput(session, "sel_mod",
                        label = "Select a module to display",
                        choices = c("No modules available" = ""),
                        selected = "")
      output$network_visual <- renderPlot({
        plot.new()
        text(0.5, 0.5, "No modules are available to visualize.")
      })
      updateTabsetPanel(session, "tabSelect", selected = "panel2")
      done(session)
      return(NULL)
    }
    
    module_choices <- seq_along(grp_idx)
    updateSelectInput(session, "sel_mod",
                      label = "Select a module to display",
                      choices = module_choices,
                      selected = module_choices[1])
    
    output$network_visual <- renderPlot({
      # # A temp file to save the output.
      # # This file will be removed later by renderImage
      # outfile <- tempfile(fileext = '.png')
      # # Generate the PNG
      # png(outfile, width = 800, height = 600)
      
      selected_module <- as.numeric(input$sel_mod)
      validate(need(is.finite(selected_module) && selected_module >= 1 &&
                      selected_module <= length(grp_idx),
                    "Select a module to display."))
      idx <- grp_idx[selected_module]
      pos_node <- which(fit_lou$membership==idx) ## position of nodes in this module
      validate(need(length(pos_node) > 0, "The selected module has no nodes to display."))
      adj_single <- dag_sub[pos_node, pos_node, drop = FALSE]
      prob_single <- mat_post_sub[pos_node, pos_node, drop = FALSE]
      nc_positions <- if (p2 > 0) seq_len(p2) else integer(0)
      a1 <- length(intersect(pos_node, nc_positions)) ## number of ncRNAs in this module
      
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
    
    
    updateTabsetPanel(session, "tabSelect", selected = "panel3")
    done(session)
  })
  
  # end
  
}
