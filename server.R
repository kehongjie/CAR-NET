# Setting the maximum file upload limit to 100 MB
options(shiny.maxRequestSize=100*1024^2)
shinyServer(function(input, output, session) {
  installed=installed.packages()[,"Package"]

  callModule(setting_server, "setting")
  ImProxy <- callModule(preproc_server, "preproc")
  callModule(saved_data_server, "saved_data", ImProxy)
  #callModule(analysis_server, "analysis")
  callModule(pathway_server, "pathway")
  callModule(indPathway_server, "indPathway")
  #callModule(combnPathway_server, "combnPathway")
  
#  if (TOOLSET.de %in% installed)
#    callModule(meta_de_server, "meta_de")
#  if (TOOLSET.clust %in% installed)
#    callModule(meta_clust_server, "meta_clust")
#  if (TOOLSET.path %in% installed)	
#    callModule(meta_path_server, "meta_path")
#  if (TOOLSET.dcn %in% installed)
#    callModule(meta_dcn_server, "meta_dcn")

  # callModule(meta_pca_server,"meta_pca")
})
