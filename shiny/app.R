library(shiny)
library(ggplot2)
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    (theme_bw(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text.x  = element_text(angle=45, hjust=1), 
                axis.line.x = element_line(colour="black"),
                axis.line.y = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                legend.direction = "vertical",
                legend.box = "vetical",
                legend.key.size= unit(0.5, "cm"),
                #legend.margin = unit(0, "cm"),
                legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
}
?theme_foundation

ui <- fluidPage(
    titlePanel("Gene expression patterns in offspring of obese dams"),
    sidebarLayout(
        sidebarPanel(
            selectizeInput('gene_selection',"Gene of interest",choices=NULL,multiple=TRUE,width="100%")
        ),
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

server <- function(input, output, session) {
    df_results <- read.csv("df_results.csv.gz")
    updateSelectizeInput(session, 'gene_selection', choices=unique(df_results$symbol), server=TRUE)
    
    output$distPlot <- renderPlot({
            ggplot(data=df_results[df_results$symbol %in% input$gene_selection,],
                   aes(x=GEOSET, y=logFC, ymin=CI.L, ymax=CI.R))+
            geom_pointrange()+
            geom_hline(yintercept=0, lty=2)+
            coord_flip() +  # flip coordinates (puts labels on y axis)
            xlab("GEOSET") + ylab("Log Fold Change (95% CI)") +
            theme_Publication()+
            ggtitle(input$gene_selection)
        
    })
}

shinyApp(ui = ui, server = server)