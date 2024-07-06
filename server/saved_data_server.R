saved_data_server <- function(input, output, session) {
  
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
  
  observeEvent(input$ACS_ADS, {
    wait(session, "Generating Bayesian Network, may take a while")
    path_old <- getwd()
    try({
      path_old <- getwd()
      filePath <- input$ncRNA_file
      print(filePath)
      
      x <- read.csv("./data/TCGA_KIRP_early_lncRNA_top5per.csv", header=T)
      rownames(x) <- x[,1]
      x <- x[,-1]
      X <- t(x)
      y <- read.csv("./data/TCGA_KIRP_early_gene_top5per.csv", header=T)
      rownames(y) <- y[,1]
      y <- y[,-1]
      Y <- t(y)

      fit <- bn.main(X=X, Y=Y, alpha=1e-5)
      
      save(fit, file="bn.RData")
      print("Bayesian Network saved as 'bn.RData'.")
      
      # output  global ACS/ADS
      output$globalACS_ADSTable <- renderImage({
        outfile <- tempfile(fileext = '.png')
        png(outfile, width = 800, height = 600)
        
        plot(graph_from_adjacency_matrix(fit$adj), mode = "directed")
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
