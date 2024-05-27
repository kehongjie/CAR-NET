rm(list=ls())
# Download R Shiny from github, extract it to a folder and rename to "NGRN_Rshiny". 
# Set the working directory to be the path to the folder "NGRN_Rshiny".
setwd("/")
shiny::runApp('NGRN_Rshiny', port=9987, launch.browser=T)


