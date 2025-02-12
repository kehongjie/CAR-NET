preproc_server <- function(input, output, session) {
  
  ns <- NS("preproc")
  
  ##########################
  # Reactive Values        #
  ##########################
  DB   <- reactiveValues(names=DB.ls(db))
  names_rna <- reactiveVal()  # final list of lncRNA lncipedia symbols
  names_gene <- reactiveVal() # final list of gene HGNC symbols
  ncRNA_name <- reactiveVal() # naming format of original lncRNA data, 
  # 0 = LNCipedia, 1 = ENSG , 2 = HSALN, 3 = HGNC
  RNA_type <- reactiveVal()   # type of original RNA data (micro or long non-coding)
  gene_name <- reactiveVal()  # naming format of original gene data
  lnc_mat <- reactiveVal()
  rna_mat <- reactiveVal()
  
  ##########################
  # Observers              #
  ##########################
  
  # watch for tab change, get the newest list of all data
  observeEvent(input$tabChange, {DB$names <- DB.ls(db)}, 
               label="tab change")
  
  # mean filter function
  internal_filter = function(M, fence){
    df = M
    avg_df = apply(df, 1, mean)
    ind = avg_df < fence
    #return(list(df[,!ind, drop = FALSE], ind))
    return(list(df[!ind, ], ind))
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
  
  
  ## lncRNA submit button -----
  observeEvent(input$submit_btn, {
    filePath <- input$ncRNA_file
    fileText <- read.csv(filePath$datapath, check.names = F, row.names = 1)
    df <- fileText
    mean_threshold <- as.numeric(input$mean)
    variance_threshold <- as.numeric(input$variance)
    zero_threshold <- input$zero
    log2 <- input$log2_transformation
    val2 <- input$RNA_type
    val <- input$ncRNA_name
    plat <- input$platform
    
    if(plat == 2 ){
      df_f <- internal_filter_zero(df,input$zero)
      #df_f <- as.data.frame(internal_filter(df, mean_threshold)[1], check.names = FALSE)
      df_f <- internal_filter(df, mean_threshold)[[1]]
      df_f <- internal_filter_var(df,variance_threshold) 
    }else if(plat == 1){
      #df_f <- as.data.frame(internal_filter(df, mean_threshold)[1], check.names = FALSE)
      df_f <- internal_filter(df, mean_threshold)[[1]]
      df_f <- internal_filter_var(df,variance_threshold) 
    }
    
    if(log2 == 2){
      df_f <- log2(df_f+1)
    }
    
    df <- t(df_f)
    genes <- rownames(df_f)
    if(val2 == 2){ # lncRNA
      ENSG_ids <- c()
      if(val == 0){ # LNCipedia
        ENSG_ids <- colnames(df)
      }else if (val == 1) { # ENSG
        ref <- read.csv("./data/reference/ENSG_lnci.csv")
        for (gene in genes) {
            if (gene %in% ref$ENSG) {
                ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
            } else {
                ENSG_ids <- c(ENSG_ids, gene)
            }
         }
      }else if (val == 2){
            ref <- read.csv("./data/reference/HSAL_lnci.csv")
            for (gene in genes) {
                if (gene %in% ref$HSAL) {
                    ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
                } else {
                    ENSG_ids <- c(ENSG_ids, gene)
          }
        }
      }else if (val == 3){
        ref <- read.csv("./data/reference/HGNC_lnci.csv")
        if (gene %in% ref$HGNC) {
          ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
        }else {
          ENSG_ids <- c(ENSG_ids, gene)
        }
      }
      colnames(df) <- ENSG_ids
      names_rna(ENSG_ids)
    }
    vals$filtered_ncRNA_data <- t(df)
    output$ncRNA_file <- DT::renderDataTable({
      vals$filtered_ncRNA_data
    })
  }, label="Filter lncRNA data.")
  
  
# watch for filter for gene expression------------
  
  observeEvent(input$submit_btn_gene, {
    #browser() 
    filePath <- input$gene_file
    fileText <- read.csv(filePath$datapath, check.names = FALSE,row.names = 1)
    df <- fileText
    mean_threshold_gene <- as.numeric(input$gene_mean)
    variance_threshold_gene <- as.numeric(input$gene_variance)
    zero_threshold_gene <- input$gene_zero
    val <- input$name_gene
    plat_g <- input$gene_platform
    log2_g <- input$gene_log2_transformation
    ref <- read.csv("./data/reference/ENSG_HGNC.csv")
    genes <- rownames(df)
    if(val == 1 ){
      ENSG_ids <- c()
      for (gene in genes) {
        if (gene %in% ref$Ensembl_gene_ID) {
          ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
        } else {
          ENSG_ids <- c(ENSG_ids, gene)
        }
      }
      rownames(df) <- ENSG_ids
    }
    
    if(plat_g == 2 ){
      df_f <- internal_filter_zero(df,zero_threshold_gene)
      df_f <- internal_filter(df, mean_threshold_gene)[[1]]
      df_f <- internal_filter_var(df,variance_threshold_gene) 
    }else if(plat_g == 1){
      df_f <- internal_filter(df, mean_threshold_gene)[[1]]
      df_f <- internal_filter_var(df,variance_threshold_gene) 
    }
    if(log2_g == 2){
      df_f <- log2(df_f+1)
    }
    vals$filtered_gene_data <- df_f
    output$gene_file <- DT::renderDataTable({
      vals$filtered_gene_data
    })
  },label="Filter gene expression data")
  
 

  
  ######################################
  # FILE UPLOAD-----
  
  # watch for ncRNA expression files upload-----
  # observeEvent(input$ncRNA_file, {
  #   if (!is.null(input$ncRNA_file)) {
  #     output$ncRNA_file <- DT::renderDataTable({
  #       filePath <- input$ncRNA_file
  #       fileText <- read.csv(filePath$datapath, check.names = FALSE)
  #       df <- t(fileText)
  #       colnames(df) <- df[1,]
  #       df <- df[-1,]
  #       df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
  #       #df <- as.data.frame(internal_filter(df, input$mean)[1], check.names = FALSE)
  #       genes <- colnames(df)
  #       names_rna(genes)
  #       val <- isolate(input$name_ncRNA)
  #       ncRNA_name(val)
  #       val2 <- isolate(input$RNA_type)
  #       RNA_type(val2)
  #       
  #       if (val2 == 2) {
  #         ENSG_ids <- c()
  #         if (val == 0) {
  #           ENSG_ids <- colnames(df)
  #         } else if (val == 1) {
  #           ref <- read.csv("./data/reference/ENSG_lnci.csv")
  #           for (gene in genes) {
  #             if (gene %in% ref$ENSG) {
  #               ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
  #             } else {
  #               ENSG_ids <- c(ENSG_ids, gene)
  #             }
  #           }
  #         } else if (val == 2) {
  #           ref <- read.csv("./data/reference/HSAL_lnci.csv")
  #           for (gene in genes) {
  #             if (gene %in% ref$HSAL) {
  #               ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
  #             } else {
  #               ENSG_ids <- c(ENSG_ids, gene)
  #             }
  #           }
  #         } else if (val == 3) {
  #           ref <- read.csv("./data/reference/HGNC_lnci.csv")
  #           for (gene in genes) {
  #             if (gene %in% ref$HGNC) {
  #               ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HGNC == gene)[1]])
  #             } else {
  #               ENSG_ids <- c(ENSG_ids, gene)
  #             }
  #           }
  #         }
  #         colnames(df) <- ENSG_ids
  #         names_rna(ENSG_ids)
  #       }
  #       t(df)
  #     })
  #   } else {
  #     return(NULL)
  #   }
  # }, label="ncRNA file upload")
  
  # watch for gene expression file upload------
  # observeEvent(input$gene_file, {
  #   if (!is.null(input$gene_file)) {
  #     output$gene_file <- DT::renderDataTable({
  #       filePath <- input$gene_file
  #       fileText <- read.csv(filePath$datapath)
  #       df <- t(fileText)
  #       colnames(df) <- df[1,]
  #       df <- df[-1,]
  #       df <- data.frame(apply(df, MARGIN = c(1,2), FUN = function(x) as.numeric(as.character(x))))
  #       #df <- as.data.frame(internal_filter(df, input$gene_mean)[1], check.names = FALSE)
  #       ref <- read.csv("./data/reference/ENSG_HGNC.csv")
  #       genes <- colnames(df)
  #       gene_name(isolate(input$name_gene))
  #       names_gene(genes)
  #       
  #       if (isolate(input$name_gene) == 1) {
  #         ENSG_ids <- c()
  #         for (gene in genes) {
  #           if (gene %in% ref$Ensembl_gene_ID) {
  #             ENSG_ids <- c(ENSG_ids, ref$Approved_symbol[which(ref$Ensembl_gene_ID == gene)])
  #           } else {
  #             ENSG_ids <- c(ENSG_ids, gene)
  #           }
  #         }
  #         colnames(df) <- ENSG_ids
  #         names_gene(ENSG_ids)
  #       }
  #       
  #       t(df)
  #     })
  #   } else {
  #     return(NULL)
  #   }
  # }, label="gene file upload")

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
