rm(list=ls())
#Download R Shiny from github and extract it to a folder "CAR-NET-main". 
#Set the working directory to be the path to the folder "CAR-NET-main".
setwd("~/Documents/GitHub/")
shiny::runApp('CAR-NET-main', port=9987, launch.browser=T)
