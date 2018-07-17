To run the shinyApp you should have R and the shiny package installed. Please download the files on your computer in a new directory and in R type: 

library(shiny)
library(rfigshare)

runApp("path_to_directory_containing_RData_folder")

This app downloads data from figshare and therefore you should have a figshare to download the needed files. 

Alternatively, you can use RStudio to visualise the shinyApp: Please download the files onto your computer in a new directory and open the ui.R file with RStudio. Then click "Run App" which is in the top right corner of the ui.R script for the application to run.
