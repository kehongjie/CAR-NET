preproc_server <- function(input, output, session) {

  ns <- NS("preproc")
  ##########################
  # Reactive Values        #
  ##########################
  DB   <- reactiveValues(names=DB.ls(db))
  STUDY <- reactiveValues(action="", 
    update=0, ori=NULL, preview=NULL,  
    expr=NULL, clinical=NULL, 
    studyName=NULL, species=NULL, genes=NULL,
    DE_p=NULL, DE_lfc=NULL, MCMC=NULL)
  SUMMARY <- reactiveValues(studySummary=data.frame(NULL), 
    DESummary=data.frame(NULL), DESummaryPM=data.frame(NULL))

  ##########################
  # Validation             #
  ##########################
  validate.study <- function(study) {
    if(is.null(STUDY$MCMC))  {
      stop(MSG.datasetInput.noinput)
    }
    studyName <- STUDY$studyName
    if(is.null(studyName) || studyName == "") {
      stop(MSG.study.noname)
    }
    if(studyName %in% DB$names) stop(MSG.study.duplicate(studyName))  
  }

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
  # NAME SELECTION
  
  # watch for name selection for ncRNA
  name_ncRNA <- 0
  
  observeEvent(input$name_ncRNA, {
    if (input$name_ncRNA == 1) {
      name_ncRNA <- 1
    } else if (input$name_ncRNA == 2) {
      name_ncRNA <- 2
    } else if (input$name_ncRNA == 3) {
      name_ncRNA <- 3
    }
  }, label="name for gene")

  # watch for name selection for gene
  name_gene <- 0
  
  observeEvent(input$name_gene, {
    if (input$name_gene == 1) {
      name_gene <- 1
    }
  }, label="name for gene")
  
  ######################################
  # FILTER SELECTION
  
  # watch for filter selection for ncRNA
  observeEvent(input$cutoff_ncRNA, {
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- as.data.frame(internal_filter(fileText[,-1], input$cutoff_ncRNA)[1], check.names = FALSE)
        
        genes <- colnames(df)
        
        if (name_ncRNA == 1) {
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
        } else if (name_ncRNA == 2) {
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
        } else if (name_ncRNA == 3) {
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
        
        if (name_gene == 1) {
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
        
        if (name_ncRNA == 1) {
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
        } else if (name_ncRNA == 2) {
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
        } else if (name_ncRNA == 3) {
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
        
        if (name_gene == 1) {
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
        # } else if (input$name_gene == 2) {
        #   ENSG_ids <- c()
        #   for (gene in genes) {
        #     if (gene %in% ref$Approved_symbol) {
        #       ENSG_ids <- c(ENSG_ids, ref$HGNC_ID[which(ref$Approved_symbol == gene)])
        #     } else {
        #       ENSG_ids <- c(ENSG_ids, "N/A")
        #     }
        #   }
        #   colnames(df) <- ENSG_ids
        # }
        
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
