# This script creates a local dataframe with the yeast proteome from Uniprot

library(stringr)
library(dplyr)

download.file("http://www.uniprot.org/docs/yeast.txt", "./yeast.txt")

#Here I was trying to read the header from the txt file but it's prob not worth it
#column_names <- (read.fwf("yeast.txt", widths = c(75,20,11,12,13,10), header = F, skip = 55, n=1, strip.white = T))

column_names <- c("Gene designations", "OLN", "Swiss-Prot AC", "Entry", "SGD X-ref", "Size", "3D", "CH")
column_classes <- c(rep("character",5), "numeric", rep("character",2))

yeast <- read.fwf("yeast.txt", widths = c(75,20,11,12,12,5,4,2), header = F, col.names = column_names, colClasses = column_classes, skip = 58, n=6726, strip.white = T)

#The dataset has observations (rows) with multiple gene names separated by a ;
#This line extracts rows with multiple Gene designation names in a new dataframe (NB: grepl returns a logical vector T when there is a ;)
yeast_mult<- yeast[grepl(";", yeast$Gene.designations),]

#make a vector that counts how many ; there are on each row
rep_vect <- str_count(yeast_mult$Gene.designations,";")

# split multiple Gene Designations into single vector elements
vect_names <- unlist(strsplit(yeast_mult$Gene.designations, ";"))
# and remove blanks from the elements
vect_names <- gsub(" ", "", vect_names)

#copy each rows as many times as the number of different gene names (taken from the rep_vect + 1)
yeast_expand <- yeast_mult[rep(row.names(yeast_mult), rep_vect+1), ] # +1 because the number of names is always one more than the number of ; 

#rename the Gene.designations column with the vector of single names
yeast_expand$Gene.designations <- vect_names

#drop rows with multiple names from the original "yeast" df using dplyr
yeast_single <- filter(yeast, !grepl(";",Gene.designations)) # filter keeps rows that do NOT have a ; in first column

#merge the two df and arrange rows by gene name
yeast_full <- full_join(yeast_single, yeast_expand) %>% arrange(Gene.designations) %>% tbl_df()
