setting_ui <- function(id, label = "global settings") {
  ns <- NS(id)

  tabPanel("Overview", value = id,
    tags$div(class = "homepage-wrapper",
      tags$section(class = "homepage-hero",
        h1("CAR-NET", class = "homepage-title"),
        p(
          "Construction and Analysis of noncoding RNA gene regulatory NETwork",
          class = "homepage-subtitle"
        ),
        tags$div(class = "homepage-actions",
          a(
            "CAR-NET R Shiny Tutorial",
            href = "tutorial/CAR-NET_tutorial.html",
            target = "_blank"
          ),
          a(
            "GitHub repository",
            href = "https://github.com/kehongjie/CAR-NET",
            target = "_blank"
          )
        )
      ),
      tags$section(class = "homepage-section",
        tags$div(class = "homepage-feature-grid",
          tags$article(class = "homepage-feature-card",
            h3("Upload data"),
            p("Choose the TCGA_KIRP example or upload matched ncRNA and gene expression CSV files.")
          ),
          tags$article(class = "homepage-feature-card",
            h3("Build network"),
            p("Run CAR-NET to infer ncRNA-gene regulatory structure from processed expression data.")
          ),
          tags$article(class = "homepage-feature-card",
            h3("Interpret results"),
            p("Review network modules, pathway findings, and visual summaries for biological interpretation.")
          )
        )
      ),
      tags$section(class = "homepage-section homepage-readable",
        h2("Quick Start"),
        tags$ol(class = "homepage-steps",
          tags$li(
            strong("Open Data Uploading and Preprocessing."),
            " Start from the upload tab in the main navigation."
          ),
          tags$li(
            strong("Select TCGA_KIRP example."),
            " The included example loads matched ncRNA and gene expression data."
          ),
          tags$li(
            strong("Review previews and preprocessing defaults."),
            " Confirm the data tables and keep the example-ready settings."
          ),
          tags$li(
            strong("Open Analysis and click Run CAR-NET."),
            " Build the example regulatory network from the prepared inputs."
          )
        )
      ),
      tags$section(class = "homepage-section homepage-readable",
        h2("Introduction"),
        p(
          "CAR-NET is an analytical software tool with an R Shiny-based graphical user interface for constructing noncoding RNA-gene regulatory networks from expression data."
        ),
        p(
          "The software implements an efficient two-stage hybrid Bayesian network structure learning approach followed by module detection and biology-driven downstream analysis."
        ),
        p(
          "CAR-NET is developed and maintained by ",
          a(
            "Dr. Tianzhou (Charles) Ma's group",
            href = "https://matianzhou.github.io",
            target = "_blank"
          ),
          " (Department of Epidemiology and Biostatistics, University of Maryland College Park)."
        ),
        p(
          "We recommend R >= 4.0.0, Rcpp >= 1.0.0, Shiny >= 1.0.0, and Python 3.11 with pandas and Pillow/PIL to run CAR-NET."
        )
      ),
      tags$section(class = "homepage-section homepage-workflow",
        h2("CAR-NET Workflow"),
        img(
          src = "flowchart.png",
          class = "homepage-flowchart",
          alt = "CAR-NET workflow diagram"
        )
      )
    )
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
