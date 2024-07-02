<<<<<<< HEAD
saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Bayesian Network Generation"),
               h3("Stage 1: Screening"),
               # numericInput(ns("permNumGlobal"),
               #              label = "Number of permutations",
               #              value = 100),

               actionButton(ns('ACS_ADS'), 'Run Stage 1',
                            class="btn-success",icon = icon("play")),

               hr(),
               h3("Stage 2: Score-based MCMC"),
               actionButton(ns("plotGlobalMDS"),
                            'Run Stage 2',
                            class="btn-success",
                            icon = icon("play")),
               tags$hr(),
               
               h2("Module Detection and Downstream Analysis"),
               actionButton(ns("plotGlobalMDS"),
                            'Divide into modules',
                            class="btn-success",
                            icon = icon("play")),
               selectInput(ns("measure"), label = "Select a module",
                           choices = list("Module 1" = "Fmeasure",
                                          "Module 2" = "youden",
                                          "Module 3" = "geo.mean"))

             ),
             mainPanel(
               tabsetPanel(tabPanel("Table of genome-wide c-scores and d-scores",
                                    h3("Genome-wide c-scores & d-scores"),
                                    textOutput(ns("Global_ACS_ADS_note")),
                                    DT::dataTableOutput(ns("globalACS_ADSTable"))),
                           tabPanel("Multidimensional scaling (MDS) map",
                                    h3("MDS map of all studies based on genome-wide c-scores"),
                                    plotOutput(ns("globalMdsFig")))
               )
             )
           )
  )
}
=======
saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Bayesian Network Generation"),
               h3("Stage 1: Screening"),
               # numericInput(ns("permNumGlobal"),
               #              label = "Number of permutations",
               #              value = 100),

               actionButton(ns('ACS_ADS'), 'Run Stage 1',
                            class="btn-success",icon = icon("play")),

               hr(),
               h3("Stage 2: Score-based MCMC"),
               actionButton(ns("plotGlobalMDS"),
                            'Run Stage 2',
                            class="btn-success",
                            icon = icon("play")),
               tags$hr(),
               
               h2("Module Detection and Downstream Analysis"),
               actionButton(ns("plotGlobalMDS"),
                            'Divide into modules',
                            class="btn-success",
                            icon = icon("play")),
               selectInput(ns("measure"), label = "Select a module",
                           choices = list("Module 1" = "Fmeasure",
                                          "Module 2" = "youden",
                                          "Module 3" = "geo.mean"))

             ),
             mainPanel(
               tabsetPanel(tabPanel("Table of genome-wide c-scores and d-scores",
                                    h3("Genome-wide c-scores & d-scores"),
                                    textOutput(ns("Global_ACS_ADS_note")),
                                    DT::dataTableOutput(ns("globalACS_ADSTable"))),
                           tabPanel("Multidimensional scaling (MDS) map",
                                    h3("MDS map of all studies based on genome-wide c-scores"),
                                    plotOutput(ns("globalMdsFig")))
               )
             )
           )
  )
}
>>>>>>> 62ec91a8086b7b7b8bff08ac6025a36d3fb03e96
