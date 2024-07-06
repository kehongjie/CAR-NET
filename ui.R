installed <- installed.packages()[,"Package"]

shinyUI(
  navbarPage("NGRN", id="nav",
    header=tagList(
      tags$div(id="working-dir",
        tagList(
          tags$p("Working Directory", class="header-label"),
          tags$p(textOutput("setting-working.dir"))
        )
      )
    ),
    # tab for global settings
    setting_ui("setting"),
    # tab for preprocessing
    preproc_ui("preproc"),
    # tab for manipulating saved data
    saved_data_ui("saved_data"),
    #analysis_ui("analysis"),

    # tab for toolsets
    #toolsets,
    tags$div(
      tags$div(id="loading",
        tags$div(id="loadingcontent",
          tags$p(id="loadingspinner", "loading......")
        )
      )
    ),
    # Including css and javascripts in head section
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "css/messenger.css"),
      tags$link(rel = "stylesheet", type = "text/css", href = "css/messenger-theme-future.css"),
      tags$link(rel = "stylesheet", type = "text/css", href = "css/styles.css"),
      tags$script(src="js/spin.min.js"),
      tags$script(src="js/messenger.min.js"),
      tags$script(src="js/message-handler.js")
    )
  )
)
