preproc_server <- function(input, output, session) {

  ns <- NS("preproc")
  ##########################
  # Reactive Values        #
  ##########################
  DB   <- reactiveValues(names=DB.ls(db))

  ##########################
  # Observers              #
  ##########################
  # watch for tab change, get the newest list of all data
  observeEvent(input$tabChange, {DB$names <- DB.ls(db)}, 
               label="tab change")
  
  internal_filter = function(M, fence){
    df = M
    avg_df = apply(df, 2, mean)
    ind = avg_df < fence
    return(list(df[,!ind, drop = FALSE], ind))
  }
  
  ######################################
  # FILTER SELECTION
  ncRNA_name <- reactiveVal()
  RNA_type <- reactiveVal()
  gene_name <- reactiveVal()
  
  # watch for filter selection for ncRNA
  observeEvent(input$cutoff_ncRNA, {
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- as.data.frame(internal_filter(fileText[,-1], input$cutoff_ncRNA)[1], check.names = FALSE)
        
        genes <- colnames(df)
        
        val <- ncRNA_name()
        val2 <- RNA_type()
        
        if (val2 == 2) {
          if (val == 1) {
            ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/ENSG_lnci.csv")
            ENSG_ids <- c()
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, "N/A")
              }
            }
            colnames(df) <- ENSG_ids
          } else if (val == 2) {
            ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/HSAL_lnci.csv")
            ENSG_ids <- c()
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, "N/A")
              }
            }
            colnames(df) <- ENSG_ids
          } else if (val == 3) {
            ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/HGNC_lnci.csv")
            ENSG_ids <- c()
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, "N/A")
              }
            }
            colnames(df) <- ENSG_ids
          }
        }
        
        df
      })
    } else {
      return(NULL)
    }
  }, label="Filter Cutoff for ncRNA")
  
  # watch for filter selection for gene
  observeEvent(input$cutoff_gene, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- as.data.frame(internal_filter(fileText[,-1], input$cutoff_gene)[1], check.names = FALSE)
        
        ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        
        if (gene_name() == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, "N/A")
            }
          }
          colnames(df) <- ENSG_ids
        }
        
        df
      })
    } else {
      return(NULL)
    }
  }, label="Filter Cutoff for ncRNA")
  
  ######################################
  # FILE UPLOAD
  
  # watch for ncRNA expression files upload
  observeEvent(input$ncRNA_file, {
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        
        df <- as.data.frame(internal_filter(fileText[,-1], input$cutoff_ncRNA)[1], check.names = FALSE)
        
        genes <- colnames(df)
        val <- isolate(input$name_ncRNA)
        ncRNA_name(val)
        val2 <- isolate(input$RNA_type)
        RNA_type(val2)
        
        if (val2 == 2) {
          if (val == 1) {
            ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/ENSG_lnci.csv")
            ENSG_ids <- c()
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, "N/A")
              }
            }
            colnames(df) <- ENSG_ids
          } else if (val == 2) {
            ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/HSAL_lnci.csv")
            ENSG_ids <- c()
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, "N/A")
              }
            }
            colnames(df) <- ENSG_ids
          } else if (val == 3) {
            ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/HGNC_lnci.csv")
            ENSG_ids <- c()
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, "N/A")
              }
            }
            colnames(df) <- ENSG_ids
          }
        }
        
        df
      })
    } else {
      return(NULL)
    }
  }, label="ncRNA file upload")
  
  # watch for gene expression file upload
  observeEvent(input$gene_file, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath)
        
        df <- as.data.frame(internal_filter(fileText[,-1], input$cutoff_gene)[1], check.names = FALSE)
        
        ref <- read.csv("C:/Users/xiaod/Downloads/Rshiny/NGRN_Rshiny/data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        gene_name(isolate(input$name_gene))
        
        if (isolate(input$name_gene) == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, "N/A")
            }
          }
          colnames(df) <- ENSG_ids
        }
        
        df
      })
    } else {
      return(NULL)
    }
  }, label="gene file upload")

  # watch for clinical file upload
  observeEvent(input$clinical_file, {
    if (!is.null(input$clinical_file)) {
      output$clinical_file <- DT::renderDataTable({
        filePath <- input$clinical_file
        fileText <- read.csv(filePath$datapath)
        fileText
      })
    } else {
      return(NULL)
    }
  }, label="clinical file upload")
}
