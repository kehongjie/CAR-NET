# CAR-NET
CAR-NET is a RShiny-based application with graphical user interface (GUI) for inferring non-coding RNA regulatory network from transcriptomic data and curated database. This application took input of both ncRNA and gene expression data generated from either bulk or single-cell RNA-seq and incorporated a list of downloaded curated databases for users to choose based on the type of ncRNA and the condition/disease of interest. It includes a preparatory step to preprocess the expression data following the standard pipeline and major analytical steps to construct the ncRNA regulatory network, identify the differential network, detect network modules, and perform pathway analysis to facilitate the biological interpretation of the network findings. In addition, it provides visualization of the network/modules and downloadable graphical and tabular outputs. 

![Alt text](./flowchart.png)

## Below is the instruction for installing and running NGRN/CAMO
#### Install the Shiny software
1. Install the CAMO R package and dependency packages following the instruction at [https://github.com/CAMO-R/Rpackage](https://github.com/CAMO-R/Rpackage).
2. Download the CAMO Shiny project at [https://github.com/CAMO-R/Rshiny](https://github.com/CAMO-R/Rshiny) by clicking on "code > Download ZIP" and extract to a local folder renamed as "Rshiny".

#### Start the Shiny software
1. Open "RunShiny.R" file in R console.
2. Set the working directory of R to the directory "path\_to\_Rshiny/" which contains the Shiny project folder "Rshiny" saved above.
3. Click on the "Run App" button.
Note that the installation progress of R packages may take up to a few minutes. Please check the progress in R console and may need to select whether to update all/some/none packages. After all packages has been installed, the CAMO Shiny app will automatically open in your default browser.


