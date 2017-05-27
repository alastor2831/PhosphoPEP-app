
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

# Retrive table of all knock-out strains in PPEP --------------------------

scraped_ko_strains <- htmlTreeParse("http://www.sbeams.org/devDC/sbeams/cgi/Glycopeptide//browse_kinases.cgi", useInternalNodes = T)
raw_ko_strains <- xpathSApply(scraped_ko_strains, "//tr[@bgcolor='#E0E0E0']/* | //tr[@bgcolor='#F1F1F1']/*", xmlValue)
index_vect <- rep(1:3, len = (length(raw_ko_strains)))
ko_strains <- data.frame(raw_ko_strains[index_vect==1], raw_ko_strains[index_vect==2], raw_ko_strains[index_vect==3])
names(ko_strains) <- c("KO Strain", "N° regulated proteins", "N° peptides")

# Scraping Phosphopep knock-outs database -----------------------------------------------------

ppep_url <- "http://www.sbeams.org/devDC/sbeams/cgi/Glycopeptide//kinase_details.cgi?kinase=" #backbone of the url(s) to be scraped
xpath_expr <- "//table[@id='regulates']/tr/*" # xpath expression that positions inside the table with the data
#urls <- character()

get_ko_data <- function(k, url = ppep_url, df = ko_strains){
    # creates a vector with the systematic name of the knock-out strains to be scraped 
    sysnames <- sys_name(k)
    ko_df_list <- vector("list") # data from each gene knock-out strain will be stored in a df inside a list
    # loop over the systematic names vector
    for (i in 1:length(sysnames)) { 
      if (sysnames[[i]] %in% df$KO == TRUE){ # check if a knock-out strain for that acc-num is present in ppep
        scraped_ko_data <- htmlTreeParse(paste0(ppep_url, sysnames[[i]]), useInternalNodes = T)
        raw_ko_data <- xpathSApply(scraped_ko_data, xpath_expr, xmlValue)
        index_vect_2 <- c(rep(0,4), rep(1:4, length.out = (length(raw_ko_data)-4))) #an index vector to build the df
        ko_df_list[[i]] <- data.frame(raw_ko_data[index_vect_2==1], raw_ko_data[index_vect_2==2], raw_ko_data[index_vect_2==3], raw_ko_data[index_vect_2==4])
        names(ko_df_list[[i]]) <- raw_ko_data[1:4]
        names(ko_df_list)[i] <- k[i]
        ko_df_list[[i]] %>% select(-KO_Strain) # drops the column KO_strain since it has the same value for each observation
        
      } else {
        ko_df_list[[i]] <- paste("The knock-out strain", k[i], "is not present in PhosphoPep")
      }
    }
    return(ko_df_list)
}
