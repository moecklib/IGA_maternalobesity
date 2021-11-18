library(shiny)

ui <- fluidPage(
    titlePanel("Obesity"),
    sidebarLayout(
        sidebarPanel(
            selectizeInput('gene_selection',NULL,choices=NULL,multiple=TRUE,width="100%")
        ),
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

server <- function(input, output, session) {
    x <- read.csv("df_results.csv.gz")
    updateSelectizeInput(session, 'gene_selection', choices=unique(x$symbol), server=TRUE)
    
    output$distPlot <- renderPlot({
        plot(0,0,type="n");text(0,0,input$gene_selection)
    })
}

shinyApp(ui = ui, server = server)
