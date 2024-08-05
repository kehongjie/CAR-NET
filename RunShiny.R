rm(list=ls())
#Download R Shiny from github and extract it to a folder "NGRN_Rshiny". 
#Set the working directory to be the path to the folder "NGRN_Rshiny".
setwd("~/Documents/GitHub/")
shiny::runApp('NGRN', port=9987, launch.browser=T)
