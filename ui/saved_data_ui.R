saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Analysis"),
               h3("Part A:"),
               h3("Screening and Score-based MCMC"),

               selectInput(ns("condition"), label = "Choose a Constraint",
                           choices = c(1,2,3)),
               selectInput(ns("database"), label = "Select a Prior Database",
                           choices = c(1,2,3)),
               selectInput(ns("tuning"), label = "Advanced Option Tuning",
                           choices = c(1,2,3)),
               actionButton(ns('ACS_ADS'), 'Generate Bayesian Network',
                            class="btn-success",icon = icon("play")),
               hr(),
               
               h3("Part B:"),
               h3("Module Detection and Downstream Analysis"),
               actionButton(ns("plotGlobalMDS"),
                            'B1: Divide into modules',
                            class="btn-success",
                            icon = icon("play")),
               h3(),
               actionButton(ns("pathway"),
                            'B2: Pathway Analysis',
                            class="btn-success",
                            icon = icon("play")),
               h3(),
               actionButton(ns("differential"),
                            'B3: Differential Modules',
                            class="btn-success",
                            icon = icon("play")),
               hr(),
               
               selectInput(ns("measure"), label = "Select a Module",
                           choices = c(1))

             ),
             mainPanel(
               tabsetPanel(id = "tabSelect",
                 tabPanel(value = "panel1",
                     h3("Inital Bayesian Network"),
                     verbatimTextOutput(ns("text")),
                     plotOutput(ns("globalACS_ADSTable"))),
                 tabPanel(value = "panel2",
                     h3("Individual Modules"),
                     plotOutput(ns("globalMdsFig")))
               )
             )
           )
  )
}
