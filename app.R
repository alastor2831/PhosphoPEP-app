knockout_list <- readRDS("./data/knockout_list.Rds")
oln_lookup <-  readRDS("./data/OLN_lookup.Rds")

library(shiny)
library(magrittr)
library(rvest)
library(DT)
library(shinythemes)
library(httr)
library(dplyr)
library(colorRamps)
library(shinyWidgets)

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

ui <- fluidPage(theme = shinytheme("slate"),
   
   h1("PhosphoPEP Database"),
   
   sidebarLayout(
      sidebarPanel(width = 3,
        h4("Select a kinase knock-out strain and press Search"),
        selectInput(inputId = "gene_name",
                      label = NULL,
                      choices = knockout_list$gene 
                     ),
        actionBttn(inputId = "button1",
                   label = "Search",
                   style = "simple",
                   color = "success",
                   size = "md",
                   block = F),
        div(style="padding-top:15px"),
        tags$hr(style="border-top: dashed 1px;"),
        h4("Get more info on the selected protein from Unirpot"),
        actionBttn(inputId = "uniprot_btn",
                   label = "Get Uniprot Data",
                   style = "simple",
                   color = "success",
                   size = "sm",
                   block = F),
        div(style="padding-bottom:10px"),
        htmlOutput(outputId = "prot_info"),
        div(style="padding-bottom:10px"),
        tableOutput(outputId = "residues_info"),
        div(class='credits',
          p("The phospho-proteomics data displayed in this website are from Phosphopep v2.0, a project by the Aebersold group at the Swiss Federal Institute of Technology (ETH) in collaboration with the Functional Genomics Center ( FGCZ ) in Zurich, Switzerland, and the Institute for Systems Biology (ISB) in Seattle, Washington USA"),
          p('Reference:'),
          HTML("<p style='font-style:italic'> Bernd Bodenmiller, Johan Malmstrom, Bertran Gerrits, David Campbell, Henry Lam, Alexander Schmidt, Oliver Rinner, Lukas N Mueller, Paul T Shannon, Patrick G Pedrioli, Christian Panse, Hoo-Keun Lee, Ralph Schlapbach, Ruedi Aebersold   Molecular Systems Biology. <b>2007</b> Oct; 3 (139)."),
          tags$b(p(HTML('<a href="http://www.nature.com/msb/journal/v3/n1/full/msb4100182.html">PhosphoPep—a phosphoproteome resource for systems biology research in Drosophila Kc167 cells</a>')),
          HTML('<div style="position:fixed;padding-bottom:15px;bottom:0;color: white;font-family: Courier;font-size: 14px;font-style: italic;text-decoration-line: underline;text-decoration-style: dashed;text-decoration-color: #28b78d;">
            <p>this is a Shiny App coded by Nick ~circa end of 2018</p></div>')
        )
      )
    ),
  
    mainPanel(
      htmlOutput(outputId = "tableTitle"),
      tabPanel("tab1", DT::DTOutput(outputId = "ppep_table"))
      )
    )
)



server <- function(input, output) {
  
  title <- eventReactive(input$button1, {
    h3_text <- HTML(paste("<h3>Currently displaying data for <i>S. cerevisiae</i><b>", input$gene_name,
                        "</b>knock-out strain from the PhosphoPEP database</h3>"))
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
      mutate(GeneNames = paste(gene, collapse = "; "),
             FoldChange = as.numeric(FoldChange)) %>%
      filter(row_number() == 1) %>% 
      ungroup() %>% 
      select(-gene, -KO_Strain) %>%
      rename("UniprotAC" = starts_with("Swi")) %>%
      mutate(UniprotAC = paste0("<a href='https://www.uniprot.org/uniprot/", UniprotAC,
                          "' target='_blank'>", UniprotAC, "</a>")) %>%
      select(GeneNames, OLN, FoldChange, Peptide, starts_with("Uni")) %>% 
      separate(GeneNames, paste0(rep("Gene", max(str_count(.$GeneNames, ";"))+1),
                                 1:max(str_count(.$GeneNames, ";")+1)),
               sep = ";", fill = "right") %>% 
      mutate_at(vars(starts_with("Gene")), funs(replace(., is.na(.), "")))
  }) 

  more_info <- eventReactive(input$uniprot_btn, {
    get_uniprot_data(as.character(oln_lookup[oln_lookup$gene == input$gene_name, 3]))
  })
  
  output$residues_info <- renderTable({
    more_info()[[3]]
  })
  
  output$prot_info <- renderText({
     paste("<p>", unlist(more_info()[1:2]), "</p>")
  })
  
  output$tableTitle <- renderText({
    title()
  })
  
  output$ppep_table <- renderDT({
   
    brks <- c(seq(min(df()$FoldChange), to = 0, by = 1.5),
              seq(0, max(df()$FoldChange), by = 1.5))
    
    length_palette <- length(brks) + 1
    negative.length <- round(abs(range(df()$FoldChange)[1]) / 
                               diff(range(df()$FoldChange)) * 
                               length_palette)
    
    positive.length <- length_palette - negative.length
    
    cols <- c(colorRampPalette(c("orangered2", "white"))(negative.length),
              colorRampPalette(c("white", "lawngreen"))(positive.length))
    
    DT::datatable(df(), style = "bootstrap", options = list(paging = FALSE,
                                                            selection = "none"),
                  escape = FALSE) %>%
      formatStyle("FoldChange",
                  backgroundColor = styleInterval(cuts = brks,
                                                  values = cols),
                  fontWeight = "bold",
                  color = "black")
    
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

