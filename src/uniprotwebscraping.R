# define helper function to get more data about the selected protein from Uniprot 

get_uniprot_data <- function(uniprot_ac){
  
  uniprot_data <- GET(paste0("https://www.uniprot.org/uniprot/", uniprot_ac, ".xml")) %>%
    content(as = "raw", content = "text/xml") %>%
    read_html
  
  protein_info <- uniprot_data %>% 
    html_nodes(xpath = '//recommendedname/* |
               //name[@type="primary"] | //comment[@type="function"]/text |
               //comment[@type="interaction"]/text') %>%
               {list(name = html_name(.), text = html_text(.))} 
  
  protein_info <- map2_chr(protein_info[[1]], protein_info[[2]], paste, sep = ": ")
  
  biogrid_link <- uniprot_data %>% 
    html_nodes(xpath = '//dbreference[@type="BioGrid"] | //dbreference[@type="BioGrid"]/*') %>%
    xml_attrs() %>%
    unlist %>%
    c(paste0(.[1], ": ", "<a href='https://thebiogrid.org/",
             .["id"],
             "' target='_blank'>ID ",
             .["id"], "</a>",
             "  (# interactions: ",
             .["value"], ")")) %>% 
    .[-(1:4)] %>% 
    as.vector()
  
  
  residues_vector <- uniprot_data %>%
    html_nodes(xpath = '//feature[@type="modified residue"] |
               //feature[@type="modified residue"]/location/*') %>%
    xml_attrs() %>%
    unlist()
  
  mod_aa <- tibble(
    type = residues_vector[seq(1, to = length(residues_vector), by = 4)],
    description = residues_vector[seq(2, to = length(residues_vector), by = 4)],
    position = residues_vector[seq(4, to = length(residues_vector), by = 4)]
  )
  
  return(list(protein_info, biogrid_link, mod_aa))
  
}