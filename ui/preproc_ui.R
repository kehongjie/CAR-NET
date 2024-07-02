<<<<<<< HEAD
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
        
        h3("ncRNA Expression Data"),
        p("All lncRNA names/IDs will be converted to LNCipedia format."),
        radioButtons(ns("RNA_type"), label = "Select the type of ncRNA data.",
                     choices = list("miRNA" = 1, 
                                    "lncRNA" = 2),
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
        p("All gene names/IDs will be converted to the HGNC approved symbol."),
        radioButtons(ns("name_gene"), label = "Select the naming format of your gene data.",
                     choices = list("HGNC Symbol" = 0, 
                                    "ENSG ID" = 1),
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
=======
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
        
        h3("ncRNA Expression Data"),
        p("All lncRNA names/IDs will be converted to LNCipedia format."),
        radioButtons(ns("RNA_type"), label = "Select the type of ncRNA data.",
                     choices = list("miRNA" = 1, 
                                    "lncRNA" = 2),
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
        p("All gene names/IDs will be converted to the HGNC approved symbol."),
        radioButtons(ns("name_gene"), label = "Select the naming format of your gene data.",
                     choices = list("HGNC Symbol" = 0, 
                                    "ENSG ID" = 1),
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
>>>>>>> 62ec91a8086b7b7b8bff08ac6025a36d3fb03e96
