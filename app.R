knockout_list <- readRDS("./data/knockout_list.Rds")
oln_lookup <-  readRDS("./data/OLN_lookup.Rds")
library(shiny)
library(tidyverse)
library(rvest)
library(DT)
library(shinythemes)

# Define UI for application that draws a histogram
ui <- fluidPage(#theme = shinytheme("slate"),
  shinythemes::themeSelector(),
   # Application title
   titlePanel("PhosphoPEP Database"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
          selectInput(inputId = "gene_name", label = "Choose a knock-out strain from the PhosphoPEP Database",
                     choices = knockout_list$gene 
                     ),
          actionButton(inputId = "button1", label = "Search")
      ),
      
      mainPanel(
        DT::DTOutput(outputId = "ppep_table")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #knock_out_table <- eventReactive(input$button1, {
    
    
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
    df()
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

