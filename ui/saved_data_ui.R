saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Analysis"),
               h3("Part A: running CAR-NET"),
               
               selectInput(ns("tune_method"), label="How to choose the tunning parameter alpha",
                           choices=list("auto-tune (will take some time!!)"=1, 
                                        "choose by myself"=2)),
               conditionalPanel(
                 condition = "input['saved_data-tune_method'] == 2",
                 numericInput(ns("alpha"), label ="Specify alpha (carefully!!)",
                              min = 0, max = 1, 
                              value = 1e-5)
               ),
               
               # selectInput(ns("condition"), label = "Choose a Constraint",
               #             choices = c(1,2,3)),
               selectInput(ns("database"), 
                           label="Select external databases to be used as prior (multiple):",
                           choices = list("miRCancer"=1,
                                          "miRTarBase"=2,
                                          "EVLncRNA"=3, 
                                          "LncRNA2Target"=4,
                                          "LncTarD"=5),
                           multiple=TRUE),
               # selectInput(ns("tuning"), label = "Advanced Option Tuning",
               #             choices = c(1,2,3)),
               # selectInput(ns("used_sample"), label="Select the samples to be used:",
               #             choices=list("all samples"=1, 
               #                          "partial samples (for differential regulation)"=2)),
               checkboxInput(ns("diff_reg"), "Want to look into differential regulation?", 
                             FALSE),
               conditionalPanel(
                 condition = "input['saved_data-diff_reg'] == true",
                 textInput(ns("cond_var"), 
                           label="Enter the name of the biological condition variable that you want to work on:", 
                           placeholder = "Enter variable name..."),
                 textInput(ns("baseline"), label="Enter the reference level:",
                           placeholder="Enter reference level name...")
               ),
               actionButton(ns('run_algo'), 'Run CAR-NET',
                            class="btn-success",icon = icon("play")),
               selectInput(ns("pathway_select"),
                           label="Select pathways to work on (multiple):",
                           choices = list("GO"=1,
                                          "KEGG"=2,
                                          "Reactome"=3,
                                          "Biocarta"=4),
                           multiple=TRUE),
               actionButton(ns("pathway_all"),
                            'B2: Do pathway analysis for all genes in the network',
                            class="btn-success",
                            icon = icon("play")),
               hr(),
               
               
               h3("Part B:"),
               h3("Module Detection"),
               numericInput(ns("mod_size"), label ="Minimum module size:",
                            min = 1, max = 100, 
                            value = 10, step=1),
               actionButton(ns("butt_mod"),
                            'B1: Detect module',
                            class="btn-success",
                            icon = icon("play")),
               # h3(),
               selectInput(ns("sel_mod"), label = "Select a module to display (in descending module size):",
                           choices = c(1)),
               # selectInput(ns("pathway"), 
               #             label="Select pathways to work on (multiple):",
               #             choices = list("GO"=1,
               #                            "KEGG"=2,
               #                            "Reactome"=3, 
               #                            "Biocarta"=4),
               #             multiple=TRUE),
               actionButton(ns("pathway_module"),
                            'B2: Do pathway analysis for genes in this module',
                            class="btn-success",
                            icon = icon("play")),
               h3(),
               # actionButton(ns("differential"),
               #              'B3: Differential Modules',
               #              class="btn-success",
               #              icon = icon("play")),
               hr(),
               
               
               h3("Part C: Differential Regulation"),
               textInput(ns("diff_level"), 
                         label="Enter another level of biological condition as opposed to the reference level:",
                         placeholder="Enter the level name..."),
               actionButton(ns("diff_run"),
                            'Update the network',
                            class="btn-success",
                            icon = icon("play")),
             ),
             
             
             mainPanel(
               tabsetPanel(id = "tabSelect",
                           # tabPanel(value = "panel1",
                           #          h3("Inital Bayesian Network"),
                           #          verbatimTextOutput(ns("text")),
                           #          plotOutput(ns("globalACS_ADSTable"))),
                           # tabPanel(value = "panel2",
                           #          h3("Individual Modules"),
                           #          plotOutput(ns("globalMdsFig")))
                           tabPanel(value = "panel1",
                                    h3("Overview of the full network"),
                                    verbatimTextOutput(ns("text")),
                                    tableOutput(ns("table_pair")),
                                    plotOutput(ns("path_plot_all"))),
                           tabPanel(value = "panel2",
                                    h3("Modulization"),
                                    plotOutput(ns("heatmap"), width = "1000px", height = "800px"),
                                    tableOutput(ns("table_node"))),
                           tabPanel(value = "panel3",
                                    h3("Modules visualization"),
                                    plotOutput(ns("network_visual"))),
                           tabPanel(value = "panel4",
                                    h3("Differential regulation"),
                                    verbatimTextOutput(ns("text4")))
                                    
               )
             )
             
             
           )
  )
}
