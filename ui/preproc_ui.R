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
        # textInput(ns("species"), label="Species name:", value = "", width = NULL,
        #           placeholder = NULL),
        # radioButtons(ns("species"), label = "Chose species",
        #         choices = list("human" = "human", "mouse" = "mouse"),
        #         selected = "human"),
        #### chose input data type
        h3("ncRNA Expression Data"),
        radioButtons(ns("name_ncRNA"), label = "Name for ncRNA",
                     choices = list("Approved Symbol" = 0,
                                    "ENSG" = 1, 
                                    "HSALN" = 2,
                                    "HGNC" = 3),
                     selected = 0),
        numericInput(ns("cutoff_ncRNA"), label = "Filter Cutoff for ncRNA",
                     value = 0,
                     min = 0,
                     max = 1,
                     step = 0.1
        ),
        p("Filter for columns whose average value is greater than this threshold."),
        
        fileInput(ns("ncRNA_file"), 'Upload ncRNA expression data file (.csv)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', 
                           '.csv')
        ),
        
        h3("Gene Expression Data"),
        radioButtons(ns("name_gene"), label = "Name for Gene",
                     choices = list("Approved Symbol" = 0,
                                    "ENSG" = 1, 
                                    "HGNC" = 2),
                     selected = 0),
        numericInput(ns("cutoff_gene"), label = "Filter Cutoff for Gene",
                     value = 0,
                     min = 0,
                     max = 1,
                     step = 0.1
        ),
        p("Filter for columns whose average value is greater than this threshold."),
        
        fileInput(ns("gene_file"), 'Upload gene expression data file (.csv)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', 
                           '.csv')
        ),
        
        # h2("Preprocess data"),
        # br(),
        # actionButton(ns('preprocSingleStudy'), 'Preprocess single study', class="btn-success"),
        # tags$hr(),

        ##########################
        # Save and Metadata      #
        ##########################
        # h2("Save single study"),
        # textInput(ns("studyName"), "Study name (Please do not use '_' in study name):", value = ""),
        # actionButton(ns('saveStudy'), 'Save', icon=icon("save"), class="btn-success")
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
