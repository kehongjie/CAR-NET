installed <- installed.packages()[,"Package"]
#enabled <- c()
#for (toolset in TOOLSET.all) {
#  if (toolset %in% installed) {
#    enabled <- c(enabled, toolset)
#  }
#}

#toolsets <- c("Toolsets")
#if (TOOLSET.de %in% enabled)
#  toolsets <- c(toolsets, list(meta_de_ui("meta_de")))
#if (TOOLSET.clust %in% enabled)
#  toolsets <- c(toolsets, list(meta_clust_ui("meta_clust")))
#if (TOOLSET.path %in% enabled)
#  toolsets <- c(toolsets, list(meta_path_ui("meta_path")))
#if (TOOLSET.dcn %in% enabled)
#  toolsets <- c(toolsets, list(meta_dcn_ui("meta_dcn")))

#if (length(toolsets) > 1) {
#  toolsets <- do.call(navbarMenu, toolsets)
#} else {
#  toolsets <- tags$div()
#}

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
    # tab for preprosessing
    preproc_ui("preproc"),
    # tab for manipulating saved data
    saved_data_ui("saved_data"),
    #analysis_ui("analysis"),
    # navbarMenu("Downstream Analysis",
    #            global_ui("global"),
    #            pathway_ui("pathway"),
    #            indPathway_ui("indPathway")
    # ),
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
