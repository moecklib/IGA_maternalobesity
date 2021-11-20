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

ui <- fluidPage(
    
    fluidRow(
        column(8,
               h2(strong("Impact of maternal obesity on gene expression")),
               h5(em("Beat Moeckli, Vaihere Delaune, Julien Prados, 
             Matthieu Tihy, Andrea Peloso, Graziano Oldani, Thomas Delmi, 
             Florence Slits, Quentin Gex, Stephanie Lacotte, Christian Toso")),
               style='margin-bottom:30px;border:1px solid; padding: 10px;'
        ),
        column(4,
               img(src = "logo_unige.jpg", height="100"))),
    
    
    fluidRow(
        column(4,
               #Determines the input gene, with possible multiple selections
               selectizeInput("gene_selection","Gene of interest",
                              choices=NULL,multiple=TRUE,width="100%"),
               
               #Select only datasets depeding on conditions
               selectizeInput(inputId = "GEOSET_type", label = strong("GEOSET type"),
                              choices = NULL, multiple=TRUE,
                              selected = "all"),
               
               checkboxGroupInput(inputId = "GEOSET", label = strong("GEOSET"),
                                  choices=NULL)
        ),
        
        column(8,
               h4("Datasets according to selection", align="center"),
               plotOutput("forrestPlot")
        )
    ),
    
    fluidRow(
        column(4,
               h4("Volcano Plot", align="center"),
               plotOutput("volcanoPlot")
        ),
        
        column(4,
               h3("Plot 3", align="center")
        ),
        
        column(4,
               h3("Plot 4", align="center")
        )
    ),
    
    fluidRow(column(12,
                    DT::dataTableOutput("table"),
                    style='margin-top:30px;border:1px solid; padding: 10px;'
    )),
    
    fluidRow(style="background-color:#f7d0e3",
             column(12,offset=0,
                    div(style="align:left",tags$small(a(icon("home",lib="glyphicon"),"Transplantation and Hepatology lab, University of Geneva, Switzerland",href="https://www.unige.ch/medecine/chiru/en/research-groups/905toso/"))),
                    div(style="align:left",tags$small(a(icon("envelope",lib="glyphicon"),"Beat Moeckli",href="mailto:beat.moeckli@etu.unige.ch"))),
                    div(style="align:left",tags$small("website designed by:",a("Beat Moeckli",href="mailto:julien.prados@unige.ch"),", ", a(href="https://www.unige.ch/medecine/bioinformatics/","Bioinformatics Support Platform, University of Geneva")))
             )
    ),
    
)



server <- function(input, output, session) {
    df_results <- read.csv("df_results.csv.gz")
    
    
    updateSelectizeInput(session, 'gene_selection', choices=unique(df_results$symbol), server=TRUE)
    
    updateCheckboxGroupInput(session, 'GEOSET', choices=unique(df_results$GEOSET))
    
    output$forrestPlot <- renderPlot({
        ggplot(data=df_results[df_results$symbol %in% input$gene_selection&
                                   df_results$GEOSET %in% input$GEOSET,],
               aes(x=GEOSET, y=logFC, ymin=CI.L, ymax=CI.R))+
            geom_pointrange()+
            geom_hline(yintercept=0, lty=2)+
            coord_flip() +  # flip coordinates (puts labels on y axis)
            xlab("GEOSET") + ylab("Log Fold Change (95% CI)") +
            theme_Publication()+
            ggtitle(input$gene_selection)
        
    })
    
    output$volcanoPlot<-renderPlot({
        #Define Unige color
        color_unige<-c("#CF0063")
        
        #Filter the 8 values from GSE46359 with outlier p-values
        volc_df<-df_results[!-log10(df_results$P.Value)>20,]
        
        #Select the gene to highlight
        df_select_genes<-volc_df[volc_df$symbol %in% input$gene_selection,]
        
        ggplot(volc_df, aes(x=logFC,y=-log10(P.Value)))+
            geom_point(alpha=0.3, size=1)+
            geom_point(data=df_select_genes, aes(x=logFC,y=-log10(P.Value)),
                       color=color_unige)+
            theme_Publication()
    })
    
    output$table <- DT::renderDataTable(DT::datatable({
        data <- df_results[df_results$symbol %in% input$gene_selection,]
        data
    }))
}

shinyApp(ui = ui, server = server)