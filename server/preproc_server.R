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
  tcga_ncRNA_path <- file.path("data", "TCGA_KIRP_early_lncRNA_top5per.csv")
  tcga_gene_path <- file.path("data", "TCGA_KIRP_early_gene_top5per.csv")
  tcga_clinical_path <- file.path("data", "TCGA_KIRP_early_clinical.csv")
  file_info_from_path <- function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }
    path <- normalizePath(path, mustWork = TRUE)
    info <- file.info(path)
    data.frame(
      name = basename(path),
      size = info$size,
      type = "text/csv",
      datapath = path,
      stringsAsFactors = FALSE
    )
  }
  file_datapath <- function(file_info) {
    file_info$datapath[1]
  }
  value_or_default <- function(value, default) {
    if (is.null(value)) default else value
  }
  selected_data_source <- function() {
    value_or_default(input$data_source, "tcga_kirp")
  }
  resolve_file_info <- function(uploaded_file, example_path) {
    if (selected_data_source() == "tcga_kirp") {
      file_info_from_path(example_path)
    } else if (!is.null(uploaded_file) && !is.null(uploaded_file$datapath)) {
      uploaded_file
    } else {
      NULL
    }
  }
  initial_ncRNA_file <- file_info_from_path(tcga_ncRNA_path)
  initial_gene_file <- file_info_from_path(tcga_gene_path)
  initial_clinical_file <- file_info_from_path(tcga_clinical_path)
  initial_clinical_data <- if (!is.null(initial_clinical_file)) {
    read.csv(file_datapath(initial_clinical_file), check.names = FALSE)
  } else {
    data.frame(
      Status = "No bundled TCGA_KIRP clinical CSV is available; clinical data are optional for this example.",
      stringsAsFactors = FALSE
    )
  }
  vals <- reactiveValues(
    file1 = initial_ncRNA_file,
    file2 = initial_gene_file,
    filtered_ncRNA_data = NULL,
    filtered_gene_data = NULL,
    clinical_file = initial_clinical_file,
    clinical_data = initial_clinical_data
  )
  
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

  process_ncRNA_file <- function(file_info, mean_threshold, variance_threshold,
                                 zero_threshold, log2_value, rna_type_value,
                                 ncRNA_name_value, platform_value) {
    if (is.null(file_info)) {
      return(NULL)
    }
    fileText <- read.csv(file_datapath(file_info), check.names = FALSE, row.names = 1)
    df <- fileText
    
    if(platform_value == 2 ){
      df_f <- internal_filter_zero(df, zero_threshold)
      #df_f <- as.data.frame(internal_filter(df, mean_threshold)[1], check.names = FALSE)
      df_f <- internal_filter(df, mean_threshold)[[1]]
      df_f <- internal_filter_var(df, variance_threshold) 
    }else if(platform_value == 1){
      #df_f <- as.data.frame(internal_filter(df, mean_threshold)[1], check.names = FALSE)
      df_f <- internal_filter(df, mean_threshold)[[1]]
      df_f <- internal_filter_var(df, variance_threshold) 
    }
    
    if(log2_value == 2){
      df_f <- base::log2(df_f+1)
    }
    
    df <- t(df_f)
    genes <- rownames(df_f)
    if(rna_type_value == 2){ # lncRNA
      ENSG_ids <- c()
      if(ncRNA_name_value == 0){ # LNCipedia
        ENSG_ids <- colnames(df)
      }else if (ncRNA_name_value == 1) { # ENSG
        ref <- read.csv("./data/reference/ENSG_lnci.csv")
        for (gene in genes) {
          if (gene %in% ref$ENSG) {
            ENSG_ids <- c(ENSG_ids, ref$name[which(ref$ENSG == gene)[1]])
          } else {
            ENSG_ids <- c(ENSG_ids, gene)
          }
        }
      }else if (ncRNA_name_value == 2){
        ref <- read.csv("./data/reference/HSAL_lnci.csv")
        for (gene in genes) {
          if (gene %in% ref$HSAL) {
            ENSG_ids <- c(ENSG_ids, ref$name[which(ref$HSAL == gene)[1]])
          } else {
            ENSG_ids <- c(ENSG_ids, gene)
          }
        }
      }else if (ncRNA_name_value == 3){
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
    t(df)
  }

  process_gene_file <- function(file_info, mean_threshold_gene,
                                variance_threshold_gene, zero_threshold_gene,
                                name_gene_value, platform_value,
                                log2_value) {
    if (is.null(file_info)) {
      return(NULL)
    }
    fileText <- read.csv(file_datapath(file_info), check.names = FALSE, row.names = 1)
    df <- fileText
    ref <- read.csv("./data/reference/ENSG_HGNC.csv")
    genes <- rownames(df)
    if(name_gene_value == 1 ){
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
    
    if(platform_value == 2 ){
      df_f <- internal_filter_zero(df, zero_threshold_gene)
      df_f <- internal_filter(df, mean_threshold_gene)[[1]]
      df_f <- internal_filter_var(df, variance_threshold_gene) 
    }else if(platform_value == 1){
      df_f <- internal_filter(df, mean_threshold_gene)[[1]]
      df_f <- internal_filter_var(df, variance_threshold_gene) 
    }
    if(log2_value == 2){
      df_f <- base::log2(df_f+1)
    }
    df_f
  }

  refresh_ncRNA_data <- function() {
    vals$file1 <- resolve_file_info(input$ncRNA_file, tcga_ncRNA_path)
    vals$filtered_ncRNA_data <- process_ncRNA_file(
      vals$file1,
      as.numeric(value_or_default(input$mean, 1)),
      as.numeric(value_or_default(input$variance, 0.1)),
      value_or_default(input$zero, 0.1),
      value_or_default(input$log2_transformation, 1),
      value_or_default(input$RNA_type, 2),
      value_or_default(input$ncRNA_name, 1),
      value_or_default(input$platform, 1)
    )
  }

  refresh_gene_data <- function() {
    vals$file2 <- resolve_file_info(input$gene_file, tcga_gene_path)
    vals$filtered_gene_data <- process_gene_file(
      vals$file2,
      as.numeric(value_or_default(input$gene_mean, 1)),
      as.numeric(value_or_default(input$gene_variance, 0.1)),
      value_or_default(input$gene_zero, 0.1),
      value_or_default(input$name_gene, 0),
      value_or_default(input$gene_platform, 1),
      value_or_default(input$gene_log2_transformation, 1)
    )
  }

  refresh_clinical_data <- function() {
    if (selected_data_source() == "tcga_kirp") {
      vals$clinical_file <- file_info_from_path(tcga_clinical_path)
    } else if (!is.null(input$clinical_file) && !is.null(input$clinical_file$datapath)) {
      vals$clinical_file <- input$clinical_file
    } else {
      vals$clinical_file <- NULL
    }
    
    if (!is.null(vals$clinical_file)) {
      vals$clinical_data <- read.csv(file_datapath(vals$clinical_file), check.names = FALSE)
    } else if (selected_data_source() == "tcga_kirp") {
      vals$clinical_data <- data.frame(
        Status = "No bundled TCGA_KIRP clinical CSV is available; clinical data are optional for this example.",
        stringsAsFactors = FALSE
      )
    } else {
      vals$clinical_data <- data.frame(
        Status = "No clinical file uploaded. Clinical data are optional.",
        stringsAsFactors = FALSE
      )
    }
  }

  load_tcga_example_data <- function() {
    ncRNA_file <- file_info_from_path(tcga_ncRNA_path)
    gene_file <- file_info_from_path(tcga_gene_path)
    
    vals$file1 <- ncRNA_file
    vals$file2 <- gene_file
    vals$filtered_ncRNA_data <- process_ncRNA_file(
      ncRNA_file,
      mean_threshold = 1,
      variance_threshold = 0.1,
      zero_threshold = 0.1,
      log2_value = 1,
      rna_type_value = 2,
      ncRNA_name_value = 1,
      platform_value = 1
    )
    vals$filtered_gene_data <- process_gene_file(
      gene_file,
      mean_threshold_gene = 1,
      variance_threshold_gene = 0.1,
      zero_threshold_gene = 0.1,
      name_gene_value = 0,
      platform_value = 1,
      log2_value = 1
    )
  }

  load_tcga_example_data()
  
  
  ######################################
  # FILTER SELECTION-----
  
# watch for filter selection for ncRNA-------
  
  ## lncRNA submit button -----
  observeEvent(input$data_source, {
    if (selected_data_source() == "tcga_kirp") {
      updateRadioButtons(session, "RNA_type", selected = 2)
      updateRadioButtons(session, "ncRNA_name", selected = 1)
      updateRadioButtons(session, "platform", selected = 1)
      updateRadioButtons(session, "log2_transformation", selected = 1)
      updateRadioButtons(session, "name_gene", selected = 0)
      updateRadioButtons(session, "gene_platform", selected = 1)
      updateRadioButtons(session, "gene_log2_transformation", selected = 1)
      load_tcga_example_data()
    } else {
      refresh_ncRNA_data()
      refresh_gene_data()
    }
    refresh_clinical_data()
  }, ignoreNULL = FALSE, label = "Select data source.")
  
  observeEvent(input$submit_btn, {
    refresh_ncRNA_data()
  }, label="Filter lncRNA data.")
  observeEvent(input$ncRNA_file, {
    refresh_ncRNA_data()
  }, ignoreNULL = FALSE, label="Load ncRNA data.")
  
  
# watch for filter for gene expression------------
  
  observeEvent(input$submit_btn_gene, {
    refresh_gene_data()
  },label="Filter gene expression data")
  observeEvent(input$gene_file, {
    refresh_gene_data()
  }, ignoreNULL = FALSE, label="Load gene expression data")
  observeEvent(input$clinical_file, {
    refresh_clinical_data()
  }, ignoreNULL = FALSE, label="Load clinical data")
  
 

  
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

  output$ncRNA_file <- DT::renderDataTable({
    validate(need(
      !is.null(vals$filtered_ncRNA_data),
      if (selected_data_source() == "tcga_kirp") {
        "TCGA_KIRP ncRNA example file is missing."
      } else {
        "Upload an ncRNA expression CSV to preview and preprocess custom data."
      }
    ))
    vals$filtered_ncRNA_data
  })
  output$gene_file <- DT::renderDataTable({
    validate(need(
      !is.null(vals$filtered_gene_data),
      if (selected_data_source() == "tcga_kirp") {
        "TCGA_KIRP gene example file is missing."
      } else {
        "Upload a gene expression CSV to preview and preprocess custom data."
      }
    ))
    vals$filtered_gene_data
  })
  output$clinical_file <- DT::renderDataTable({
    vals$clinical_data
  })
  
  # Values to be read by saved_data_server.R
  observe({vals$cutoff_ncRNA <- input$cutoff_ncRNA})
  observe({vals$cutoff_gene <- input$cutoff_gene})
  observe({vals$names_rna <- names_rna()})
  observe({vals$names_gene <- names_gene()})

  return(vals)
}
