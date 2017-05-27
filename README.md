# PPEP-scraper
an R script to retrieve data from the PhosphoPep database (https://www.ncbi.nlm.nih.gov/pubmed/17940529)

## Why?
While using the PhosphoPep web-interface to browse through differentially expressed phosphopeptides data in yeast (S. cerevisiae), I noticed that the single knock-out strains page would only show the Systematic Name (eg. YLL021W) instead of / in addition to the standard name (eg. SPA2) for each differentially expressed phosphopeptide, forcing you to click on each one to find out the standard name of the corresponding protein.

This script written in R uses a local dataframe of the yeast proteome (retrived from Uniprot database) to find the correspondence between Systematic and Standard names. The data of selected knock-out strains are scraped from PhosphoPep website, as no APIs are available.
