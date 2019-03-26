knockout_list <- readRDS("./data/knockout_list.Rds")
oln_lookup <-  readRDS("./data/OLN_lookup.Rds")
library(shiny)
library(tidyverse)
library(rvest)
library(DT)
library(shinythemes)
library(httr)
# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("slate"),
  
   # Application title
   titlePanel("PhosphoPEP Database"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
          selectInput(inputId = "gene_name", label = "Choose a knock-out strain from the PhosphoPEP Database",
                     choices = knockout_list$gene 
                     ),
          actionButton(inputId = "button1", label = "Search"),
          textOutput(outputId = "selected_strain")
      ),
      
      mainPanel(tabsetPanel(
        tabPanel("tab1", DT::DTOutput(outputId = "ppep_table"))
      )
   )
  )
  )

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$selected_strain <- renderText({
    paste("You are looking at knock-out strain for", input$gene_name, " from PPEP DB")
  })
  
  df <- eventReactive(input$button1, {
    kinase_url <- paste0("http://www.sbeams.org/devDC/sbeams/cgi/Glycopeptide//kinase_details.cgi?kinase=",
                       knockout_list[knockout_list$gene == input$gene_name, 1])
    kinase_url %>%
      read_html() %>%
      html_nodes(xpath = '//table[@id="regulates"]') %>% 
      html_table() %>%
      .[[1]] %>% 
      as.tibble() %>% 
      rename(KO_Strain = X1, OLN = X2, FoldChange = X3, Peptide = X4) %>% 
      filter(row_number() != 1) %>% 
      inner_join(oln_lookup, by = "OLN") %>%
      group_by(OLN, Peptide) %>%
      mutate(GeneNames = paste(gene, collapse = ",")) %>%
      filter(row_number() == 1) %>% 
      ungroup() %>% 
      select(-gene)# %>% 
      #separate(GeneNames, into = paste0(rep("GeneName", 4), 1:4))
  })
  
  output$ppep_table <- renderDT({
    DT::datatable(df(), style = "bootstrap",
                  options = list(paging = FALSE, selection = "none")) %>% 
      formatStyle(
        "FoldChange", color = styleInterval(0, c("red", "green"))
      )
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)


# 
# 
# 
# raw <- GET("https://www.uniprot.org/uniprot/P12345.xml") %>% 
# 
# parsed <- content(raw, as = "raw", content = "text/xml") %>% 
#   read_html()
#   
#   parsed %>% 
#     html_nodes(xpath = '//recommendedname/* |
#                //name[@type="primary"] | //comment[@type="function"]/text |
#                //comment[@type="interaction"]/text') %>% 
#     html_text()
# 
#   parsed %>% 
#     html_nodes(xpath = '//dbreference') %>% 
#     html_table()  
#   
#                    
# residues_vector <- parsed %>%
#     html_nodes(xpath = '//feature[@type="modified residue"] |
#                //feature[@type="modified residue"]/location/*') %>% 
#     xml_attrs() %>%
#     unlist() %>%
#     tibble(
#       type = .[seq(1, to = length(.), by = 4)]),
#       description = residues_vector[seq(2, to = length(residues_vector), by = 4)],
#       position = residues_vector[seq(4, to = length(residues_vector), by = 4)]
#     )
# 
# %>% View
# 
# 
# parsed %>% 
# html_structure()