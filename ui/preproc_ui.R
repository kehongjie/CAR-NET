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
        h3("ncRNA Expression Data"),
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
        
#Filtering----------------------------------------------------------------------
        
        radioButtons(ns("platform"), label = "Select the sequencing platform of ncRNA data.",
                     choices = list("bulk RNA-seq" = 1, 
                                    "microarray" = 2,
                                    "scRNA-seq" = 3
                     ),
                     selected = 1),
# bulk RNA-seq fitlering criteria        
conditionalPanel(
  condition = "input['preproc-platform'] == 1",

  radioButtons(ns("format"), label = "Select the naming format of your lncRNA data.",
               choices = list("RPKM" = 0,
                              "TPM" = 1, 
                              "CPM" = 2),
               selected = 0),
  
  numericInput(ns("Variance"), label ="Set cutoff for variance quantile (%):",
              min = 0, max = 0, value = 0, step = 1),
  
  conditionalPanel(
    condition = "input['preproc-format'] == 0",  # RPKM selected
    numericInput(ns("cutoff_RPKM"), label = "Set filtering cutoff for RPKM:",
                 value = 0,
                 min = 0,
                 step = 1
    ),
  ),
  conditionalPanel(
    condition = "input['preproc-format'] == 1",  # TPM selected
    numericInput(ns("cutoff_TPM"), label = "Set filtering cutoff for TPM:",
                 value = 0,
                 min = 0,
                 step = 1
    ),
  ),
  conditionalPanel(
    condition = "input['preproc-format'] == 2",  # CPM selected
  #   sliderInput(ns("cutoff_CPM"), label = "Set filtering cutoff for CPM:",
  #               min = 10, max = 20, value = 10, step = 1)
  # ),
  numericInput(ns("cutoff_CPM"), label = "Set filtering cutoff for CPM:",
               value = 0,
               min = 0,
               step = 1
  ),
),
),

# microarray fitlering criteria 
conditionalPanel(
  condition = "input['preproc-platform'] == 2",
  numericInput(ns("variance"), label = "Filter variance for ncRNA",
               value = 0,
               min = 0,
               step = 1
  ),
  numericInput(ns("mean"), label = "Filter mean for ncRNA",
               value = 0,
               min = 0,
               step = 1
  ),
),

# scRNA-seq filtering criteria

conditionalPanel(
  condition = "input['preproc-platform'] == 3",
  numericInput(ns("percentage of zero count"), label = "Filter zero count for ncRNA",
               value = 0,
               min = 0,
               step = 1
  ),
  numericInput(ns("variance"), label = "Filter variance for ncRNA",
               value = 0,
               min = 0,
               step = 1
  ),
  numericInput(ns("mean"), label = "Filter mean for ncRNA",
               value = 0,
               min = 0,
               step = 1
  ),
),

# log2 transformation
        radioButtons(ns("log2_transformation"), label = "Has the data been log2-transformed?",
                     choices = list("yes" = 1, 
                                    "no" = 2
                     ),
                     selected = 1),
        
        # numericInput(ns("cutoff_ncRNA"), label = "Filter Cutoff for ncRNA",
        #              value = 0,
        #              min = 0,
        #              step = 1
        # ),
        # p("Filter for rows whose average value is greater than this threshold."),
        
        h3("Gene Expression Data"),
        p("All gene names/IDs will be converted to the HGNC approved symbol."),
        radioButtons(ns("name_gene"), label = "Select the naming format of your gene data.",
                     choices = list("HGNC Symbol" = 0, 
                                    "ENSG ID" = 1),
                     selected = 0),
        fileInput(ns("gene_file"), 'Upload gene expression data file (.csv)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', 
                           '.csv')
        ),
        numericInput(ns("cutoff_gene"), label = "Filter Cutoff for Gene",
                     value = 0,
                     min = 0,
                     step = 1
        ),
        p("Filter for rows whose average value is greater than this threshold."),
        
        h3("Clinical Data"),
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
