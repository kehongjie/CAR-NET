saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Bayesian Network Generation"),
               h3("Part A:"),
               h4("Screening and Score-based MCMC"),

               actionButton(ns('ACS_ADS'), 'Generate Bayesian Network',
                            class="btn-success",icon = icon("play")),
               hr(),
               
               h3("Part B:"),
               h4("Module Detection and Downstream Analysis"),
               
               actionButton(ns("plotGlobalMDS"),
                            'Divide into modules',
                            class="btn-success",
                            icon = icon("play")),
               hr(),
               
               selectInput(ns("measure"), label = "Please Generate the Modules First",
                           choices = c(1))

             ),
             mainPanel(
               tabsetPanel(id = "tabSelect",
                 tabPanel(value = "panel1",
                     h3("Inital Bayesian Network"),
                     plotOutput(ns("globalACS_ADSTable"))),
                 tabPanel(value = "panel2",
                     h3("Individual Modules"),
                     plotOutput(ns("globalMdsFig")))
               )
             )
           )
  )
}
