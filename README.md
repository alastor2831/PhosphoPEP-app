# PPEP-scraper
an R Shiny App to retrieve and visualize data from the PhosphoPep database (https://www.ncbi.nlm.nih.gov/pubmed/17940529) 

## Why?
While using the PhosphoPep web-interface to browse differentially expressed phosphopeptides data in yeast (*S. cerevisiae*), I noticed that the single knock-out strains page would only show the Systematic Name (eg. YLL021W) instead of/in addition to the standard name (eg. SPA2) for each differentially expressed phosphopeptide. This forced me to have to click on each peptide to find out the standard name of the corresponding protein.

Since everyone knows that less clicks is better, this app shows all the informations in one table. 

## How to use it
Launch the live version of the app from          and select the protein of interest. All phosphorilation data will be displayed in a table on the right. That's it! 

More info from Uniprot about the selected protein can be retrieved by the button in the left column.


## How does it work?

The app uses two local databases and a helper function in the src folder

#### src/localdbs.R
It downloads the entire yeast proteome database from Unirpot in a text file.
After some data wrangling and cleaning it outputs a dataframe with one row for each protein and its corresponding Systematic Name and SwissProt ID.

#### src/uniprotwebscraping.R
Since Uniprot does not provide APIs to quickly access data, this function uses webscraping to collect a brunch of information about a single protein.

#### /app.R

* **ui**

Contains the front-end code, based on a dark theme ("slated") (via the `shinyWidgets` package). It creates the two click buttons and 


* **server**

Contains two main reactive events 

1. eventReactive(input$button1)
    + web scraping of the PhosphoPEP database for the requested protein
    + data wrangling and cleaning
    + returns a data frame

2. eventReactive(input$uniprot_btn)
    + calls get_uniprot_data
    + returns a list

Contains the output functions

1. output$ppep_table
  i) diplays the interactive table on the right side panel. Uses `DT` and `colorRamp` packages to color code the `foldChange` column based on its values (green: overexpression/upregulated, red: downregulated)







This script written in R uses a local dataframe of the yeast proteome (retrived from Uniprot database) to find the correspondence between Systematic and Standard names. The data of selected knock-out strains are scraped from PhosphoPep website, as no APIs are available.
