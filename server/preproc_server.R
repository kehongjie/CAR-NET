preproc_server <- function(input, output, session) {
  
  ns <- NS("preproc")
  
  ##########################
  # Reactive Values        #
  ##########################
  DB   <- reactiveValues(names=DB.ls(db))
  names_rna <- reactiveVal()  # final list of lncRNA lncipedia symbols
  names_gene <- reactiveVal() # final list of gene HGNC symbols
  ncRNA_name <- reactiveVal() # naming format of original lncRNA data
  RNA_type <- reactiveVal()   # type of original RNA data (micro or long non-coding)
  gene_name <- reactiveVal()  # naming format of original gene data

  
  ##########################
  # Observers              #
  ##########################
  
  # watch for tab change, get the newest list of all data
  observeEvent(input$tabChange, {DB$names <- DB.ls(db)}, 
               label="tab change")
  
  # mean filter function
  internal_filter = function(M, fence){
    df = M
    avg_df = apply(df, 2, mean)
    ind = avg_df < fence
    return(list(df[,!ind, drop = FALSE], ind))
  }
  
  # variance filter function
  internal_filter_var <- function(M,fence){
    # cutoff threshold: fence
    df = M
    gene_variances <- apply(df, 1, var)
    variance_cutoff <- quantile(gene_variances, 1-fence) 
    # quantiles of variance (keep the 20% genes with the highest variance)
    df_filter <- df[gene_variances>variance_cutoff,]
    return(df_filter)
  }
  
  # zero filter funtion
  internal_filter_zero <- function(raw_counts,cell_threshold){
    min_cells <- ncol(raw_counts) * cell_threshold
    nonzero_counts <- rowSums(raw_counts > 0)
    filtered_genes <- raw_counts[nonzero_counts >= min_cells, ]
    return(filtered_genes)
  }
  
  
  ######################################
  # FILTER SELECTION-----
  
# watch for filter selection for ncRNA-------
  
  observe({
    print(input$mean)
    print(input$variance)
  })
  
  ## filter for mean ------
  observeEvent(input$mean, {
    #browser()
    #mean_text <- input$mean
    #mean_num <- suppressWarnings(as.numeric(mean_text))
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- as.data.frame(internal_filter(df, input$mean)[1], check.names = FALSE)
        genes <- colnames(df)
        names_rna(genes)
        
        val <- ncRNA_name()
        val2 <- RNA_type()
        
        if (val2 == 2) {
          ENSG_ids <- c()
          if (val == 0) {
            ENSG_ids <- colnames(df)
          } else if (val == 1) {
            ref <- read.csv("./data/reference/ENSG_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 2) {
            ref <- read.csv("./data/reference/HSAL_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 3) {
            ref <- read.csv("./data/reference/HGNC_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          }
          colnames(df) <- ENSG_ids
          names_rna(ENSG_ids)
        }
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="Filter Mean Cutoff for ncRNA")
  
  ## filter for variance ------
  
  observeEvent(input$variance, {
    #variance_text <- input$variance
    #variance_num <- suppressWarnings(as.numeric(variance_text))
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- internal_filter_var(df,input$variance)
        genes <- colnames(df)
        names_rna(genes)
        
        val <- ncRNA_name()
        val2 <- RNA_type()
        
        if (val2 == 2) {
          ENSG_ids <- c()
          if (val == 0) {
            ENSG_ids <- colnames(df)
          } else if (val == 1) {
            ref <- read.csv("./data/reference/ENSG_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 2) {
            ref <- read.csv("./data/reference/HSAL_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 3) {
            ref <- read.csv("./data/reference/HGNC_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          }
          colnames(df) <- ENSG_ids
          names_rna(ENSG_ids)
        }
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="Filter Mean Cutoff for ncRNA")
  
  ## filter for zero count scRNA-seq (ncRNA)------
  observeEvent(input$zero, {
    #  browser() 
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- internal_filter_zero(df,input$zero)
        genes <- colnames(df)
        names_rna(genes)
        
        val <- ncRNA_name()
        val2 <- RNA_type()
        
        if (val2 == 2) {
          ENSG_ids <- c()
          if (val == 0) {
            ENSG_ids <- colnames(df)
          } else if (val == 1) {
            ref <- read.csv("./data/reference/ENSG_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 2) {
            ref <- read.csv("./data/reference/HSAL_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 3) {
            ref <- read.csv("./data/reference/HGNC_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          }
          colnames(df) <- ENSG_ids
          names_rna(ENSG_ids)
        }
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="Filter zero count cutoff for ncRNA")
  
  ## Log2 or not (ncRNA)-------
  observeEvent(input$log2_transformation, {
    #browser()
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        genes <- colnames(df)
        names_rna(genes)
        
        val <- ncRNA_name()
        val2 <- RNA_type()
        
        if (val2 == 2) {
          ENSG_ids <- c()
          if (val == 0) {
            ENSG_ids <- colnames(df)
          } else if (val == 1) {
            ref <- read.csv("./data/reference/ENSG_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 2) {
            ref <- read.csv("./data/reference/HSAL_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 3) {
            ref <- read.csv("./data/reference/HGNC_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          }
          colnames(df) <- ENSG_ids
          names_rna(ENSG_ids)
        }
        t(df)
        if(input$log2_transformation==2){
          t(log2(df+1))
        }
        else{
          return(t(df))
        }
        })
    } else {
      return(NULL)
    }
  }, label="Log transformation for ncRNA")
  
  
# watch for filter for gene expression------------
  
  ## mean filter selection for gene--------
  observeEvent(input$gene_mean, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- as.data.frame(internal_filter(df, input$gene_mean)[1], check.names = FALSE)
        ref <- read.csv("./data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        names_gene(genes)
        
        if (gene_name() == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, gene)
            }
          }
          colnames(df) <- ENSG_ids
          names_gene(ENSG_ids)
        }
        
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="Filter Mean Cutoff for gene expression")
  
  
  ## variance filter for gene expression------
  observeEvent(input$gene_variance, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- internal_filter_var(df,input$gene_variance)
        ref <- read.csv("./data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        names_gene(genes)
        
        if (gene_name() == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, gene)
            }
          }
          colnames(df) <- ENSG_ids
          names_gene(ENSG_ids)
        }
        
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="Filter Variance Cutoff for gene expression")
  
  
  ## zero count filter for gene expression--------
  observeEvent(input$gene_zero, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- internal_filter_zero(df,input$gene_zero)
        ref <- read.csv("./data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        names_gene(genes)
        
        if (gene_name() == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, gene)
            }
          }
          colnames(df) <- ENSG_ids
          names_gene(ENSG_ids)
        }
        
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="Filter Zero Count Cutoff for gene expression")
  
  ## log2 or not (gene expression)----
  observeEvent(input$gene_log2_transformation, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        ref <- read.csv("./data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        names_gene(genes)
        
        if (gene_name() == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, gene)
            }
          }
          colnames(df) <- ENSG_ids
          names_gene(ENSG_ids)
        }
        t(df)
        if(input$gene_log2_transformation==2){
          t(log2(df+1))
        }
        else{
          return(t(df))
        }
        })
    } else {
      return(NULL)
    }
  }, label="Log2 transformation for gene expression")
  
  
  ######################################
  # FILE UPLOAD-----
  
  # watch for ncRNA expression files upload-----
  observeEvent(input$ncRNA_file, {
   # browser() 
    if (!is.null(input$ncRNA_file)) {
      output$ncRNA_file <- DT::renderDataTable({
        filePath <- input$ncRNA_file
        fileText <- read.csv(filePath$datapath, check.names = FALSE)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- as.data.frame(internal_filter(df, input$mean)[1], check.names = FALSE)
        genes <- colnames(df)
        names_rna(genes)
        val <- isolate(input$name_ncRNA)
        ncRNA_name(val)
        val2 <- isolate(input$RNA_type)
        RNA_type(val2)
        
        if (val2 == 2) {
          ENSG_ids <- c()
          if (val == 0) {
            ENSG_ids <- colnames(df)
          } else if (val == 1) {
            ref <- read.csv("./data/reference/ENSG_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 2) {
            ref <- read.csv("./data/reference/HSAL_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HSAL) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          } else if (val == 3) {
            ref <- read.csv("./data/reference/HGNC_lnci.csv")
            for (gene in genes) {
              if (gene %in% ref$HGNC) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
              } else {
                ENSG_ids <- c(ENSG_ids, gene)
              }
            }
          }
          colnames(df) <- ENSG_ids
          names_rna(ENSG_ids)
        }
        t(df)
      })
    } else {
      return(NULL)
    }
  }, label="ncRNA file upload")
  
  # watch for gene expression file upload------
  observeEvent(input$gene_file, {
    if (!is.null(input$gene_file)) {
      output$gene_file <- DT::renderDataTable({
        filePath <- input$gene_file
        fileText <- read.csv(filePath$datapath)
        df <- t(fileText)
        colnames(df) <- df[1,]
        df <- df[-1,]
        df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
        df <- as.data.frame(internal_filter(df, input$gene_mean)[1], check.names = FALSE)
        ref <- read.csv("./data/reference/ENSG_HGNC.csv")
        genes <- colnames(df)
        gene_name(isolate(input$name_gene))
        names_gene(genes)
        
        if (isolate(input$name_gene) == 1) {
          ENSG_ids <- c()
          for (gene in genes) {
            if (gene %in% ref$Ensembl_gene_ID) {
              ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
            } else {
              ENSG_ids <- c(ENSG_ids, gene)
            }
          }
          colnames(df) <- ENSG_ids
          names_gene(ENSG_ids)
        }
        
        t(df)
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
  
  # Values to be read by saved_data_server.R
  vals <- reactiveValues()
  observe({vals$file1 <- input$ncRNA_file})
  observe({vals$file2 <- input$gene_file})
  observe({vals$cutoff_ncRNA <- input$cutoff_ncRNA})
  observe({vals$cutoff_gene <- input$cutoff_gene})
  observe({vals$names_rna <- names_rna()})
  observe({vals$names_gene <- names_gene()})
  
  return(vals)
}
