rm(list=ls())
<<<<<<< HEAD
#Download R Shiny from github and extract it to a folder "NGRN_Rshiny". 
#Set the working directory to be the path to the folder "NGRN_Rshiny".
setwd("C:\\Users\\xiaod\\Downloads\\Rshiny")
=======
# Download R Shiny from github, extract it to a folder and rename to "NGRN_Rshiny". 
# Set the working directory to be the path to the folder "NGRN_Rshiny".
setwd("/")
>>>>>>> 62ec91a8086b7b7b8bff08ac6025a36d3fb03e96
shiny::runApp('NGRN_Rshiny', port=9987, launch.browser=T)


