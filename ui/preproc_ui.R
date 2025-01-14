preproc_ui <- function(id, label= "preprocessing data") {
  ns <- NS(id)
  tabPanel("Data Uploading and Preprocessing", value=id,
    sidebarLayout(
      sidebarPanel(
        useShinyjs(),
        ##########################
        # Upload Data
        ##########################
        h2("Upload Data"),
        p("Note: It is assumed that the uploaded data will have RNA/genes on the rows, and samples on the columns."),
# ncRNA --------
        h3("NcRNA Expression Data"),
        fileInput(ns("ncRNA_file"), 'Upload ncRNA expression data file (.csv)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', 
                           '.csv')
        ),
        p("All lncRNA names/IDs will be converted to LNCipedia format."),
        radioButtons(ns("RNA_type"), label = "Select the type of ncRNA data.",
                     choices = list("miRNA" = 1, 
                                    "lncRNA" = 2,
                                    "other" = 3,
                                    "mixed" = 4),
                     selected = 1),
        conditionalPanel(
          condition = "input['preproc-RNA_type'] == 2",
          radioButtons(ns("name_ncRNA"), label = "Select the naming format of your lncRNA data.",
                       choices = list("LNCipedia" = 0,
                                      "ENSG" = 1, 
                                      "HSALN" = 2,
                                      "HGNC" = 3),
                       selected = 0),
        ),
        
#Filtering ncRNA----------------------------------------------------------------
        
        radioButtons(ns("platform"), label = "Select the sequencing platform of ncRNA data.",
                     choices = list("bulk RNA-seq / microarray" = 1, 
                                    "scRNA-seq" = 2
                     ),
                     selected = 1),
        ## bulk RNA-seq / microarray fitlering criteria  (ncRNA) -----     
        conditionalPanel(
          condition = "input['preproc-platform'] == 1",
          numericInput(ns("variance"), label ="Keep the ? genes with the highest variance:",
                      min = 0.1, max = 1, 
                      value = 0.1,
                      step = 0.05
                      ),
          
          numericInput(ns("mean"), label = "Filter mean",
                       value = 1,
                       min = 1,
                       step = 0.5
          ),
         ),

        ## scRNA-seq filtering criteria (ncRNA)-----
        
        conditionalPanel(
          condition = "input['preproc-platform'] == 2",
          numericInput(ns("zero"), label = "Filter zero count",
                       value = 0.1,
                       min = 0,
                       max = 1,
                       step = 0.05
          ),
          numericInput(ns("variance"), label = "Keep the ? genes with the highest variance:",
                       value = 0.1,
                       min = 0,
                       max = 1,
                       step = 0.05
          ),
          numericInput(ns("mean"), label = "Filter mean",
                       value = 1,
                       min = 0,
                       step = 0.5
          ),
        ),

        ## log2 transformation------
        radioButtons(ns("log2_transformation"), label = "Has the data been log2-transformed?",
                     choices = list("yes" = 1, 
                                    "no" = 2
                     ),
                     selected = 1),

#Filtering gene expression------------------------------------------------------
        h3("Gene Expression Data"),
        p("All gene names/IDs will be converted to the HGNC approved symbol."),
        fileInput(ns("gene_file"), 'Upload gene expression data file (.csv)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', 
                           '.csv')
        ),
        radioButtons(ns("name_gene"), label = "Select the naming format of your gene data.",
                     choices = list("HGNC Symbol" = 0, 
                                    "ENSG ID" = 1),
                     selected = 0),
        
        radioButtons(ns("gene_platform"), label = "Select the sequencing platform of gene expression data.",
                     choices = list("bulk RNA-seq" = 1, 
                                    "scRNA-seq" = 2
                     ),
                     selected = 1),
        ## bulk RNA-seq fitlering criteria (gene expression)----       
        conditionalPanel(
          condition = "input['preproc-gene_platform'] == 1",
          
          numericInput(ns("gene_variance"), label ="Keep the ? genes with the highest variance:",
                       min = 0, max = 1, value = 0.1, step = 0.05),
          
          numericInput(ns("gene_mean"), label = "Filter mean",
                       value = 1,
                       min = 0,
                       step = 0.5
          ),
        ),

        ## scRNA-seq filtering criteria (gene expression)-------
        
        conditionalPanel(
          condition = "input['preproc-gene_platform'] == 2",
          numericInput(ns("gene_zero"), label = "Filter zero count",
                       value = 0.1,
                       min = 0,
                       max = 1,
                       step = 0.05
          ),
          numericInput(ns("gene_variance"), label = "Keep the ? genes with the highest variance:",
                       value = 0.1,
                       min = 0,
                       max = 1,
                       step = 0.05
          ),
          numericInput(ns("gene_mean"), label = "Filter mean",
                       value = 1,
                       min = 0,
                       step = 0.5
          ),
        ),

        ## log2 transformation (gene expression)------
      radioButtons(ns("gene_log2_transformation"), label = "Has the data been log2-transformed?",
                   choices = list("yes" = 1, 
                                  "no" = 2
                   ),
                   selected = 1),
              p("Filter for rows whose average value is greater than this threshold."),

# Clinical Data--------
        h3("Clinical Data (Optional)"),
        fileInput(ns("clinical_file"), 'Upload clinical data file (.csv)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', 
                           '.csv')
        ),
      ),

      mainPanel(
        h3("ncRNA Expression"),
        DT::dataTableOutput(ns("ncRNA_file")),
        hr(),
        h3("Gene Expression"),
        DT::dataTableOutput(ns("gene_file")),
        hr(),
        h3("Clinical Data"),
        DT::dataTableOutput(ns("clinical_file"))
      )
    )
  )
}
