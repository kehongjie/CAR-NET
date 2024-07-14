# NGRN

## Development logs (for internal communication/documentation only)
- May 24 (from Hongjie): I set up our NGRN Rshiny using CAMO as a template, so I would recommend everyone first download and run the CAMO on your own computer to get a sense of how it works. After you get CAMO running, try to download and run our current version NGRN Rshiny (I pretty much only change the welcome page as of now). Next, we can all contribute, make changes accordingly and remove unnecessary parts.
- June 19 (from Eric): I created a preliminary page for the data uploading/preprocessing with the options to change the ID/name of each ncRNA and gene based on a reference data set provided by Zhangchi, and to adjust the filter cutoff for the averages of each column. Test and reference data sets are uploaded under the data folder.
- July 14 (from Eric): Over the past 6 weeks, I've created a functional data preprocessing page whose output serves as the input to the analysis page. I implemented Zhangchi's data filter code, Hongjie's bayesian network and partition code, as well as Cindy's network visualization code. The RShiny can now generate a network visualization based on input data, as well as divide the network into individual modules. I hope that this serves as a useful template for further development.

## Below is the instruction for installing and running NGRN
#### Install the Shiny software
1. Install the CAMO R package and dependency packages following the instruction at [https://github.com/CAMO-R/Rpackage](https://github.com/CAMO-R/Rpackage).
2. Download the CAMO Shiny project at [https://github.com/CAMO-R/Rshiny](https://github.com/CAMO-R/Rshiny) by clicking on "code > Download ZIP" and extract to a local folder renamed as "Rshiny".

#### Start the Shiny software
1. Open "RunShiny.R" file in R console.
2. Set the working directory of R to the directory "path\_to\_Rshiny/" which contains the Shiny project folder "Rshiny" saved above.
3. Click on the "Run App" button.
Note that the installation progress of R packages may take up to a few minutes. Please check the progress in R console and may need to select whether to update all/some/none packages. After all packages has been installed, the CAMO Shiny app will automatically open in your default browser.


