#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("backend.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("RBP-RNA interactions"),
   p("Data supported by the binding sites of RBPs derived from CLIP-seq data", a(href="https://rnasysu.com/encori/rbpClipRNA.php?source=mRNA", "(ENCORI database)"),""),
   p("Reference transcript information obtained from", a(href="https://www.ensembl.org/info/genome/genebuild/mane.html", target="_blank", "(MANE.GRCh38.v1.0.ensembl_genomic)"), ""),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectizeInput(inputId = "rbp_list",
                       label = "RBP:",
                       choices = NULL,
                       multiple = F,
                       options = list(
                         placeholder = "Select RBP",
                         options = list(create = FALSE)),
                       selected = NULL),
        shiny::checkboxInput(inputId = "all_RBPs",
                             label = "All RBPs",
                             value = F),
        selectizeInput(inputId = "gene_list",
                       label = "Gene of interest:",
                       choices = NULL,
                       multiple = F,
                       options = list(
                         placeholder = "Select gene",
                         options = list(create = FALSE)),
                       selected = NULL),
        actionButton(inputId = "mainButton", label = "Accept")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("bindingSitesPlot", height = 800) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
         DT::DTOutput("bindingSitesTable")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
  updateSelectizeInput(session, 'rbp_list', choices = RBP_choices, server = TRUE, selected = "ADAR")
  
  
  ## UPDATE DROPDOWN GENE LIST DEPENDING ON THE RBP SELECTED --------------------------------------------
  
  observeEvent(input$rbp_list, {
    updateSelectizeInput(session, 'gene_list', choices = get_gene_list_per_RBP(input$rbp_list), server = TRUE)
  })
  

  ## PRODUCE THE PLOT -----------------------------------------------------------------------------------
  
  toListen_plot <- eventReactive(list(input$mainButton,input$gene_list), {
    visualiseCLIP(target_RBP = input$rbp_list, target_gene = input$gene_list, allRBPs = input$all_RBPs)
  })
  
  output$bindingSitesPlot <- renderPlot({
    
    toListen_plot()
    
  },width = "auto", height = "auto")
  
  
  ## PRODUCE A TABLE -----------------------------------------------------------------------------------
  
  toListen_table <- eventReactive(list(input$mainButton,input$gene_list), {
    tableCLIP(target_RBP = input$rbp_list, target_gene = input$gene_list, allRBPs = input$all_RBPs)
  })
  
  output$bindingSitesTable <- DT::renderDataTable( toListen_table() , 
                                                   server = T,
                                                   extensions = c('Buttons','RowGroup','Responsive'),
                                                   options = list(pageLength = 10,
                                                                  dom = 'Bfrtip',
                                                                  buttons = c('copy', 'csv', 'excel'),
                                                                  width = "100%",
                                                                  rownames = FALSE))

  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

