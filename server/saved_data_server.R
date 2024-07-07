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
  
  # Bayesian Network Generation
  source("./bn_main.R")
  source("./gui_filter.R")
  
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
      
      save(fit, file="bn.RData")
      print("Bayesian Network saved as 'bn.RData'.")
      
      partition_matrix_blocks <- function(adj_matrix, block_size) {
        # Number of nodes
        n <- nrow(adj_matrix)
        
        # Calculate the number of blocks along one dimension
        num_blocks <- ceiling(n / block_size)
        
        # Create a list to store the submatrices
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
      
      result <- partition_matrix_blocks(fit$adj, 10)
      
      # output  global ACS/ADS
      output$globalACS_ADSTable <- renderImage({
        outfile <- tempfile(fileext = '.png')
        png(outfile, width = 800, height = 600)
        
        plot(graph_from_adjacency_matrix(result[[1]]), mode = "directed")
        dev.off()
        
        # Return a list containing the filename
        list(src = outfile,
             contentType = 'image/png',
             width = 800,
             height = 600,
             alt = "This is alternate text")
      }, deleteFile = TRUE)
      
      sendSuccessMessage(session, "Bayesian Network generated.")
    }, session)
    
    setwd(path_old)
    done(session)
  })
  
  source("./Function_early_late.R")
  
  observeEvent(input$plotGlobalMDS, {
    wait(session, "Generating Graphs")
    
    try({
      library(utils)
      
      output$globalMdsFig <- renderImage({
        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- tempfile(fileext = '.png')
        
        # Generate the PNG
        png(outfile, width = 800, height = 600)
        p=6 # determine how many ncRNA in the matrix first # first 6 rows
        early_late_igraph(adj_early,adj_late,prob_early,prob_late,ncRNA_num=p,early_graph_name="early_edge_color_graph2.jpeg",late_graph_name="late_edge_color_graph2.jpeg")
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
    done(session)
  })
}
