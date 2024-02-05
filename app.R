#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("backend_data_preparation.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("RBP-RNA interactions"),
   p("Data supported by the binding sites of RBPs derived from CLIP-seq data (ENCORI database)."),
   p("Reference transcript information obtained from 'MANE.GRCh38.v1.0.ensembl_genomic'."),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectizeInput(inputId = "rbp_list",
                       label = "RBP:",
                       choices = NULL,
                       multiple = T,
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
  updateSelectizeInput(session, 'gene_list', choices = gene_choices, server = TRUE)

  ## PRODUCE THE PLOT -----------------------------------------------------------------------------------
  
  toListen_plot <- eventReactive(list(input$mainButton), {
    visualiseCLIP(target_RBP = input$rbp_list, target_gene = input$gene_list, allRBPs = input$all_RBPs)
  })
  
  output$bindingSitesPlot <- renderPlot({
    
    toListen_plot()
    
  },
  width = "auto",
  height = "auto")
  
  
  ## PRODUCE A TABLE -----------------------------------------------------------------------------------
  
  toListen_table <- eventReactive(list(input$mainButton), {
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

