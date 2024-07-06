saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Bayesian Network Generation"),
               h3("Stage 1/2: Screening and Score-based MCMC"),

               actionButton(ns('ACS_ADS'), 'Run Stage',
                            class="btn-success",icon = icon("play")),

               hr(),
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
              tabsetPanel(tabPanel(h3("Inital Bayesian Network"),
                               plotOutput(ns("globalACS_ADSTable"))),
                            tabPanel(h3("Directed Acyclic Graph of RNA-gene relationships"),
                               plotOutput(ns("globalMdsFig"))))
             )
           )
  )
}
