# New version -----------------------------------------------------------
# In this version the occurrences (rows) with multiple Standard Names are separated into multiple columns without row duplication
# By counting ; in the first column of the yeast df, the maximum number of multiple Standard Names (or Aliases) is 9.

library(stringr)
library(dplyr)
library(tidyr)

download.file("http://www.uniprot.org/docs/yeast.txt", "./yeast.txt")
column_names <- c("Gene designations", "OLN", "Swiss-Prot AC", "Entry", "SGD X-ref", "Size", "3D", "CH")
column_classes <- c(rep("character",5), "numeric", rep("character",2))
yeast <- read.fwf("yeast.txt", widths = c(75,20,11,12,12,5,4,2), header = F, col.names = column_names, colClasses = column_classes, skip = 58, n=6726, strip.white = T)

yeast2 <- yeast %>% separate(Gene.designations, into = c("Gene.designations", paste0(rep("Alias.", 9), 1:9)), sep = ";", fill = "right") %>% select(Gene.designations, OLN:CH, Alias.1:Alias.9) %>% tbl_df


# The sys_name function is unchanged since 
sys_name <- function(x, df = yeast_full){
    output <- df[match(x, df$Gene.designations), ]$OLN
    names(output) <- x
    return(output)
}

gene_name <- function(y, df = yeast_full){
    out <- df[match(x, df$OLN), ]$Gene.designations
    
    df_temp <- df[grep(paste(y, collapse = "|"), df$OLN, value = FALSE), c(1,2,4)] %>% arrange(Entry)
    # df_temp <- df_temp[order(df$Entry),] # another way to order the df that doesn't require dplyr
    output <- df_temp$Gene.designations
    names(output) <- df_temp$OLN
    return(output)
}
