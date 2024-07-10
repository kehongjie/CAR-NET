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
  
  
  ##########################
  # Observers              #
  ##########################
  
  adj_mat <- reactiveVal(value = NULL)
  p <- reactiveVal()
  q <- reactiveVal()
  size <- reactiveVal()
  
  # Bayesian Network Generation
  source("./bn_main.R")
  source("./gui_filter.R")
  source("./single-adj_different_size.r")
  
  observeEvent(input$ACS_ADS, {
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
      adj_mat(fit$adj)
      
      save(fit, file="bn.RData")
      print("Bayesian Network saved as 'bn.RData'.")
      
      output$text <- renderText({ paste("After Stage 1:\n",
        "Left with ", fit$s1_lvl1, " level-1 edges\n",
        "Left with ", fit$s1_lvl2, " level-2 edges\n",
        "Involving ", fit$s1_rna, " ncRNAs and ", fit$s1_gene, "genes\n",
        "After Stage 2:\n",
        "Left with ", fit$s2_lvl1, " level-1 edges\n", 
        "Left with ", fit$s2_lvl2, " level-2 edges\n",
        "Involving ", fit$s2_rna, " ncRNAs and ", fit$s2_gene, " genes", sep="") })
      
      small <- fit$adj[rowSums(fit$adj) > 0, colSums(fit$adj) > 0]
      
      # output  global ACS/ADS
      output$globalACS_ADSTable <- renderImage({
        outfile <- tempfile(fileext = '.png')
        png(outfile, width = 800, height = 600)
        
        # plot(graph_from_adjacency_matrix(result[[1]]), mode = "directed")
        p(ncol(X))
        q(ncol(Y))
        # p(108)
        # q(462)
        total <- p() + q()
        single_igraph(adj_single=fit$adj, prob_single=matrix(0.5, total, total), 
                      ncRNA_num=p(),
                      single_graph_name="./plots/test_all.png") # call your visualization function
        dev.off()
        
        # Return a list containing the file name
        list(src = outfile,
             contentType = 'image/png',
             width = 800,
             height = 600,
             alt = "This is alternate text")
      }, deleteFile = TRUE)
      
      sendSuccessMessage(session, "Bayesian Network generated.")
    }, session)
    
    observe({
      updateTabsetPanel(session, "tabSelect",
                        selected = "panel1")
    })
    done(session)
  })
  
  source("./Function_early_late.R")
  
  partition_matrix_blocks <- function(adj_matrix, block_size) {
    # Number of nodes
    n <- nrow(adj_matrix)
    
    # Calculate the number of blocks along one dimension
    num_blocks <- ceiling(n / block_size)
    
    # Create a list to store the sub matrices
    blocks <- vector("list", num_blocks^2)
    
    # Initialize the index for the list
    index <- 1
    
    # Fill the submatrices
    for (i in seq(1, n, by = block_size)) {
      for (j in seq(1, n, by = block_size)) {
        # Determine the rows and columns for the current block
        row_indices <- i:min(i + block_size - 1, n)
        col_indices <- j:min(j + block_size - 1, n)
        
        # Extract the submatrix
        submatrix <- adj_matrix[row_indices, col_indices]
        
        # Store the submatrix in the list
        blocks[[index]] <- submatrix
        index <- index + 1
      }
    }
    
    # Return the list of submatrices
    return(blocks)
  }
  
  observeEvent(input$plotGlobalMDS, {
    wait(session, "Generating Graphs")
    
    try({
      library(utils)
      
      size(235)
      result <- partition_matrix_blocks(adj_mat(), size())
      
      l <- as.list(seq(1,2*(p()+q())/size()))
      observe({
        updateSelectInput(session, "measure",
                          label = "Select a module",
                          choices = l,
                          selected = l[1])
      })
      
      output$globalMdsFig <- renderImage({
        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- tempfile(fileext = '.png')
        # Generate the PNG
        png(outfile, width = 800, height = 600)
        # p=6 # determine how many ncRNA in the matrix first # first 6 rows
        # early_late_igraph(adj_early,adj_late,prob_early,prob_late,ncRNA_num=p,early_graph_name="early_edge_color_graph2.jpeg",late_graph_name="late_edge_color_graph2.jpeg")
        # single_igraph(adj_single=result[[1]], prob_single=matrix(0.5, p()+q(), p()+q()), 
        #               ncRNA_num=p(),
        #               single_graph_name="./plots/modulex.png") # call your visualization function
        # plot(graph_from_adjacency_matrix(result[[1]]), mode = "directed")
        total <- size()
        single_igraph(adj_single=result[[as.numeric(input$measure)]],
                      prob_single=matrix(0.5, total, total), 
                      ncRNA_num=p()/2,
                      single_graph_name="./plots/test_all.png")
        dev.off()
        
        # Return a list containing the filename
        list(src = outfile,
             contentType = 'image/png',
             width = 800,
             height = 600,
             alt = "This is alternate text")
      }, deleteFile = TRUE)
      
      sendSuccessMessage(session, "Graphs are generated.")
    }, session)
    
    observe({
      updateTabsetPanel(session, "tabSelect",
                        selected = "panel2")
    })
    done(session)
  })
  
  observeEvent(input$measure, {
    if (!is.null(adj_mat())) {
      wait(session, "Generating Graphs")
      
      try({
        partitions <- partition_matrix_blocks(adj_mat(), size())
        output$globalMdsFig <- renderImage({
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          # Generate the PNG
          png(outfile, width = 800, height = 600)
          # plot(graph_from_adjacency_matrix(partitions[[as.numeric(input$measure)]]), mode = "directed")
          total <- size()
          single_igraph(adj_single=partitions[[as.numeric(input$measure)]],
                        prob_single=matrix(0.5, total, total), 
                        ncRNA_num=p()/2,
                        single_graph_name="./plots/test_all.png")
          dev.off()
    
          # Return a list containing the filename
          list(src = outfile,
               contentType = 'image/png',
               width = 800,
               height = 600,
               alt = "This is alternate text")
        }, deleteFile = TRUE)
      }, session)
      done(session)
    }
  }, label = "module selection")
}
