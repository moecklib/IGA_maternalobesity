#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Libraries and overall functions####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Things to do in the scrip:
#-Add plot which shows the expression level for each gene
#-Add a tile plot with the following selections 1)p-value, 2)logFC, 3)number of times in different datasets 4)head 50
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


#*#*#*#*#*#*#*#*#*#
#User Interface####
#*#*#*#*#*#*#*#*#*#
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
    
    #Row with forrest plot
    fluidRow(
        column(4,
               #Determines the input gene, with possible multiple selections
               selectizeInput("gene_selection","Gene of interest",
                              choices=NULL,multiple=TRUE,width="100%"),
               
               #Select only datasets depeding on conditions (sex, offspring diet, age)
               selectInput(inputId = "GEOSET_sex", label = strong("Sex offspring"),
                              choices = NULL, multiple=FALSE,
                              selected = "all"),
               selectizeInput(inputId = "GEOSET_diet", label = strong("Diet offspring"),
                              choices = NULL, multiple=FALSE,
                              selected = "all"),
               selectizeInput(inputId = "GEOSET_age", label = strong("Age offspring"),
                              choices = NULL, multiple=FALSE,
                              selected = "all")
             
        ),
        
        column(8,
               h4("Datasets according to selection", align="center"),
               plotOutput("forrestPlot")
        )
    ),
    
    
    
    fluidRow(
        column(6,
               h4("Volcano Plot", align="center"),
               plotOutput("volcanoPlot")
        ),
        
        column(6,
               h3("Plot 3", align="center")
        )
    ),
    
    #Layout of the row with the tile plot
    fluidRow(
        column(4,
               #Select gene selection criteria
               numericInput(inputId = "pValTile", label = strong("Max p-value"),
                           min = 0.0001, max=1, value=0.05),
               
               sliderInput(inputId = "logFCTile", label = strong("Min log fold change"),
                           min = 0, max=3, value=0.1, step=0.1),
               
               sliderInput(inputId = "nbTile", label = strong("Number of datasets"),
                           min = 0, max=12, value=3)
               
        ),
        
        column(8,
               h4("Genes according to selection", align="center"),
               plotOutput("tilePlot")
        )
    ),
    
    #Insertion of the datatable
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

#*#*#*#*#*#*#*#*#*#*#
#Server Functions####
#*#*#*#*#*#*#*#*#*#*#
server <- function(input, output, session) {
    #Import of dataset and selection of appropriate variables
    df_results <- read.csv("df_results.csv.gz")[,c(1,3:7, 9:11)]
    
    #Import description data of the different GEOSET
    GEOSET_descr<-local({ GEOSET_descr<-data.frame(
        GEOSET=unique(df_results$GEOSET),
        sex=c(c(rep("male",times=9), "female", "male", "male")),
        age=c("pre-natal",c(rep("adult",times=8)), "suckling", "suckling", "pre-natal"),
        offspring_diet=c(NA, "ND", c(rep("HFD",times=3)), c(rep(c("ND", "HFD"),times=2)), "ND", "ND", NA)
    )
    
    row_all<-data.frame(
        GEOSET=unique(df_results$GEOSET), 
        sex=c(rep("all",times=12)), age=c(rep("all",times=12)), offspring_diet=c(rep("all",times=12))
    )
    
    GEOSET_descr<-rbind(row_all,GEOSET_descr)
    })
    
    #update select for forrest plot
    updateSelectizeInput(session, 'gene_selection', choices=unique(df_results$symbol), server=TRUE)
    
    updateSelectInput(session, 'GEOSET_sex', choices=unique(GEOSET_descr$sex))
    
    updateSelectizeInput(session, 'GEOSET_diet', choices=unique(GEOSET_descr$offspring_diet), server=TRUE)
    
    updateSelectizeInput(session, 'GEOSET_age', choices=unique(GEOSET_descr$age), server=TRUE)
    
    #update inputs for the heatplot
    updateNumericInput(session, inputId = "pValTile")
    updateSliderInput(session, inputId = "logFCTile")
    updateSliderInput(session, inputId = "nbTile")
    
    #Create Forrest Plot
    output$forrestPlot <- renderPlot({
        
        #Select GEOSET's according to above input
        sex<-GEOSET_descr$GEOSET[GEOSET_descr$sex%in%input$GEOSET_sex]
        age<-GEOSET_descr$GEOSET[GEOSET_descr$age%in%input$GEOSET_age]
        diet<-GEOSET_descr$GEOSET[GEOSET_descr$offspring_diet%in%input$GEOSET_diet]
        
        GEO_sel<-Reduce(intersect, list(sex,age,diet))
        
        #Create plot
        ggplot(data=df_results[df_results$symbol %in% input$gene_selection&
                                   df_results$GEOSET %in% GEO_sel,],
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
    
    output$barPlot<-renderPlot({
        #Define Unige color
        color_unige<-c("#CF0063")
        
        ggplot()+
            geom_point(alpha=0.3, size=1)+
            geom_point(data=df_select_genes, aes(x=logFC,y=-log10(P.Value)),
                       color=color_unige)+
            theme_Publication()
    })
    
    #Create tile plot for display on the third row
    output$tilePlot<-renderPlot({
        ##Produce top gene list for tile plot
        #Select genes with specified criteria
        pre_tile<-df_results$symbol[abs(df_results$logFC)>input$logFCTile&
                                        df_results$P.Value<input$pValTile]
        genes_tiles<-pre_tile[table(pre_tile)>input$nbTile]
        
        #Producedf with top 40 genes
        df_genes_tile<-df_results[df_results$symbol%in%genes_tiles,]
        df_genes_tile<-df_genes_tile[order(abs(df_genes_tile$logFC), decreasing=TRUE),]
        genes_tile40<-unique(head(df_genes_tile, 45)[,2])
        df_tile40<-df_results[df_results$symbol%in%genes_tile40,]
        
        #convert high logFC (smaller or larger than 2) to 2 and select appropriate columns
        df_tile40$logFC<-ifelse(df_tile40$logFC>2, 2,df_tile40$logFC)
        df_tile40$logFC<-ifelse(-2>df_tile40$logFC, -2,df_tile40$logFC)
        df_tile40<-df_tile40[,c(1:3)]
        
        #Produce the plot to display
        ggplot(df_tile40, aes(x = GEOSET, y = reorder(symbol, logFC), fill=logFC)) +   #reorder variables according to values
            scale_fill_viridis_c()+
            labs(x="GEO Expression Set", y="Gene Symbol", fill="logFC")+
            geom_tile(colour="white",size=0.1)+
            theme_Publication()
    })
    
    output$table <- DT::renderDataTable(DT::datatable({
        data <- df_results[df_results$symbol %in% input$gene_selection,c(1:3, 6:8)]
        #data <- data[.,c()]
        data
    }))
}

shinyApp(ui = ui, server = server)