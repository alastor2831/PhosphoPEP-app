library(stringr)
library(dplyr)
library(magrittr)


download.file("http://www.uniprot.org/docs/yeast.txt", "./yeast.txt")

# Here I was trying to read the header from the txt file but it's prob not worth it
# column_names <- (read.fwf("yeast.txt", widths = c(75,20,11,12,13,10), header = F, skip = 55, n=1, strip.white = T))

column_names <- c("gene", "OLN", "Swiss-Prot AC", "Entry", "SGD X-ref", "Size", "3D", "CH")
column_classes <- c(rep("character",5), "numeric", rep("character",2))

yeast <- read.fwf("yeast.txt", widths = c(75,20,11,12,12,5,4,2),
                  header = F, col.names = column_names, colClasses = column_classes,
                  skip = 58, n=6726, strip.white = T) %>%
  as.tibble()

# save on file the Uniprot Accession Number <> OLN lookup table

yeast %>%
  separate_rows(gene, sep = '; ') %>% 
  select(1,2,3) %>% 
  saveRDS("./data/OLN_lookup.Rds") 

# save on file the data.frame of all knock-out strains in PhosphoPEP

url_knockouts <- "http://www.sbeams.org/devDC/sbeams/cgi/Glycopeptide//browse_kinases.cgi"

url_knockouts %>%
  read_html() %>% 
  html_nodes(xpath = "//table[@border=0]") %>%
  .[4] %>% 
  html_table() %>% 
  .[[1]] %>% 
  as.tibble %>% 
  rename(OLN = X1, RegulatedProteins = X2, Peptides = X3) %>% 
  filter(row_number() != 1) %>% 
  saveRDS("./data/knockout_list.Rds")
