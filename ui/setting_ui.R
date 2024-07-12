setting_ui <- function(id, label = "global settings") {
  ns <- NS(id)
  
  tabPanel("Overview", value=id,
           h2("Welcome to NGRN",align = "middle",style="primary"),
           tags$hr(),
           #HTML('<center><img src="Rshiny_flowchart.png" width="500" </center>'),
           fluidRow(
             h3("Introduction"),
             p("NGRN is an analytical software with R-Shiny based graphical user interface (GUI) tool to construct noncoding RNA-gene regulatory network from expression data."),
             p("The software implemnts our proposed efficient two-stage hybrid BN structure learning approach followed by module detection and biology-driven downstream analysis."), 
             p("Our tool is available for download on github: ",a(strong("NGRN"), href="https://github.com/kehongjie/NGRN",target="_blank"), ".", 
             "For detailed implementation of the tool, please refer to our ",a(strong("Tutorials."), href="https://github.com/metaOmicsCAMO/tutorial/
blob/master/CAMO_turtorial.pdf",target="_blank")), 
             p("NGRN is developed and maintained by ", a("Dr. Tianzhou (Charles) Ma's group ",href="https://matianzhou.github.io",target="_blank"), "(Department of Epidemiology and Biostatistics, University of Maryland College Park)."),
             p("We recommend users to use R (>=3.5.0) to implement our tool. If you are using R 3.4, 
               you may encounter errors in installing dependencies of the modules. "),
             hr(),
             h3("Instructions"),
             p("Upload your data onto the Data Uploading and Preprocessing Page, and the corresponding output will serve as the input to the analysis page."),
             style="text-indent: 20px; font-size: 16px; margin-left: 20px"),
           img(src='flowchart.png',align="middle",width="1000"),
    #        tags$hr(),
    # mainPanel(
    #   h2("Session Information"),
    #   verbatimTextOutput(ns("urlText")),
    #   h2("Saving directory:", style="display:inline"),
    #   helpIcon("working_dir_help", "During the computation, some output files or images are automatically saved to this directory."),
    #   directoryInput(ns('directory'), label='select a directory')
    # )
  )
}
