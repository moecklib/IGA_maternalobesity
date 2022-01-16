#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Libraries and overall functions####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

library(shiny)
library(ggplot2)
theme_Publication <- function(base_size=18, base_family="sans") {
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
                legend.position = "bottom",
                legend.direction = "horizontal",
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
#Define plotwidth
plotwidth<-"700px"

ui <- fluidPage(
  tags$head(includeHTML(("google-analytics.html"))),
    
    fluidRow(
        column(8,
               h1(strong("Impact of maternal obesity on gene expression")),
               h4(em("Beat Moeckli, Vaihere Delaune, Julien Prados, 
             Matthieu Tihy, Andrea Peloso, Graziano Oldani, Thomas Delmi, 
             Florence Slits, Quentin Gex, Laura Rubbia-Brandt, Nicolas Goossens, 
                     Stephanie Lacotte, Christian Toso")),
               h4("Geneva University Hospitals & Geneva University Faculty of Medicine"),
               style='margin-bottom:30px;padding: 10px;'
        ),
        
        column(4,
               img(src = "logo_unige.jpg", height="100")
               )
        ),
    
    #Row with forrest plot
    fluidRow(column(12,
        tabsetPanel(
            tabPanel(strong("Forest Plot"),
                     fluidRow(
                         column(4,
                                #Determines the input gene, with possible multiple selections
                                selectizeInput("gene_selection","Gene of interest",
                                               choices=NULL,multiple=TRUE),
                                
                                #Select only datasets depeding on conditions (sex, offspring diet, age)
                                selectInput(inputId = "GEOSET_sex", label = strong("Sex offspring"),
                                            choices = NULL, multiple=FALSE,
                                            selected = "all"),
                                selectizeInput(inputId = "GEOSET_diet", label = strong("Diet offspring"),
                                               choices = NULL, multiple=FALSE,
                                               selected = "all"),
                                selectizeInput(inputId = "GEOSET_age", label = strong("Age offspring"),
                                               choices = NULL, multiple=FALSE,
                                               selected = "all"),
                                p("Enter the gene of interest in the first field and select 
                                the characetistics of the datasets 
                                in which you would like to to display the log fold change." )
                         ),
                         
                         column(8,
                                plotOutput("forrestPlot",width = plotwidth)
                         )
                     ),
                     
                     #Insertion of level of expression plot
                     fluidRow(style="margin-top:20px",
                              column(4,
                                     h4(strong("Description")),
                                     p("The graph to the right displays the gene expression level of the 
                                     above selected gene in each individual dataset as percentile. Highly
                                       expressed genes are closer to 100%.")
                                     ),
                              
                              column(8,
                                     plotOutput("exprPlot",width = plotwidth)
                                     )
                              ),

                     
                     #Insertion of the datatable
                     fluidRow(column(4,
                                            h4(strong("Description")),
                                            p("The table to the right contains the summary result of the differential
                                              gene expression anaylsis per condition analysed (GEOSET)
                                              for the above selected gene.")
                     ),
                     
                     column(8,
                                     dataTableOutput("table"),
                                     style='margin-top:20px;padding: 10px;'
                                     )
                     )
                     
            ),
            
            #TabPanel heatmap
            tabPanel(strong("Heatmap"),
                     #Layout of the row with the tile plot
                     column(4,
                            #Select gene selection criteria
                            numericInput(inputId = "pValTile", label = strong("Max p-value"),
                                         min = 0.0001, max=1, value=0.05),
                            
                            sliderInput(inputId = "logFCTile", label = strong("Min log fold change"),
                                        min = 0, max=3, value=0.1, step=0.1),
                            
                            sliderInput(inputId = "nbTile", label = strong("Number of datasets"),
                                        min = 0, max=12, value=3),
                            
                            sliderInput(inputId = "nbGenes", label = strong("Number of genes"),
                                        min = 1, max=100, value=40),
                            
                            p("You can select genes based on the above conditions:" ),
                            p(strong("p-value:"),"Maximal p-value for the gene to be included"),
                            p(strong("logFC:"),"Minimum logarithmic fold change for the gene to be included"),
                            p(strong("Number of datasets:"),
                              "Minimum number of conditions/datasets in which the above conditions
                              need to be fulfilled for the gene to be included."),
                            p(strong("Number of genes to display:"),
                              "Number of genes that are displayed in the heatmap.")
                            
                     ),
                     
                     column(8,
                            h4("Genes according to selection", align="center"),
                            plotOutput("tilePlot", width = plotwidth, height=plotwidth)
                     ),
                     
            ),
            
            tabPanel(strong("Volcano Plot"),
                     fluidRow(
                       column(4,
                              #Determines the input gene, with possible multiple selections
                              selectizeInput("gene_selectionV","Gene of interest",
                                             choices=NULL,multiple=TRUE),
                              
                              p("Enter the gene of interest in the field above. The black diamond
                                indicates the relative position (logFC, p-value) of the selcted  genes in 
                                regards to all other genes" )
                       ),
                       
                       column(8,
                              plotOutput("volcanoPlot",width = plotwidth, height=plotwidth)
                       )
                     ),
                     
            ),
            
            #TabPanel GEOSET description
            tabPanel(strong("GEOSET description"),
                     column(4,
                            h4(strong("Description")),
                            p("The table to the right contains all informations regarding
                            the included datasets in this study with the corresponding
                            references.")
                     ),
                     
                      column(8,
                                     dataTableOutput("GEOSETtable")
                             )
            )
        )
        ),
        
        #Footer of the page
        fluidRow(
             column(12,offset=0,
                    div(style="align:left",tags$small(a(icon("home",lib="glyphicon"),"Transplantation and Hepatology lab, University of Geneva, Switzerland",
                                                        href="https://www.unige.ch/medecine/chiru/en/research-groups/905toso/"))),
                    div(style="align:left",tags$small(a(icon("envelope",lib="glyphicon"),"Beat Moeckli",href="mailto:beat.moeckli@etu.unige.ch"))),
                    div(style="align:left",tags$small("website designed by:",a("Beat Moeckli",href="mailto:beat.moeckli@etu.unige.ch"),", ", 
                                                      a(href="https://www.unige.ch/medecine/bioinformatics/","Bioinformatics Support Platform, University of Geneva"))),
                    style='margin-top:5px;padding: 20px;'
             )
    )
)
)

#*#*#*#*#*#*#*#*#*#*#
#Server Functions####
#*#*#*#*#*#*#*#*#*#*#

server <- function(input, output, session) {
    #Import of dataset and selection of appropriate variables
    #df_results <- read.csv("df_results.csv.gz")[,c(1,3:7, 9:10)];saveRDS(df_results,file="df_results.rds")
    df_results <- readRDS("df_results.rds")
    
    #Import description data table for tab of GEOSET description
    GEOSET_descr_tab<-read.csv("description_GEOSET.csv")
        
    
    # Compute expression rank
    exp_rank <- reactive({
        A <- tapply(df_results$AveExpr,list(df_results$symbol,df_results$GEOSET),mean,na.rm=TRUE)
        R <- apply(A,2,rank,ties.method="first") - 1
        R <- t(t(R) / (colSums(!is.na(A)) - 1))
        R[R>1] <- NA
        R  
    })
    
    #Import description data of the different GEOSET
    GEOSET_descr<-local({ 
        GEOSET_descr <- data.frame(
            GEOSET=unique(df_results$GEOSET),
            sex=c(c(rep("male",times=9), "female", "male", "male")),
            age=c("pre-natal",c(rep("adult",times=8)), "suckling", "suckling", "pre-natal"),
            diet_offspring=c(NA, "Normal Diet", c(rep("High Fat Diet",times=3)), 
                             c(rep(c("Normal Diet", "High Fat Diet"),times=2)), c(rep("Normal Diet",times=2)), NA)

        )
    
        row_all<-data.frame(
            GEOSET=unique(df_results$GEOSET), 
            sex=c(rep("all",times=12)), age=c(rep("all",times=12)), 
            diet_offspring=c(rep("all",times=12))
        )
        
        GEOSET_descr<-rbind(row_all,GEOSET_descr)
    })
    
    #update select for forrest plot
    updateSelectInput(session, 'GEOSET_sex', choices=unique(GEOSET_descr$sex))
    
    updateSelectizeInput(session, 'GEOSET_diet', choices=unique(GEOSET_descr$diet_offspring), 
                         server=TRUE)
    
    updateSelectizeInput(session, 'GEOSET_age', choices=unique(GEOSET_descr$age), server=TRUE)
    
    updateSelectizeInput(session, 'gene_selection', choices=unique(df_results$symbol),
                         selected = "Lcn2", server=TRUE, options=list(maxItems=1))
    
    #Gene Selection for Volcano Plot
    updateSelectizeInput(session, 'gene_selectionV', choices=unique(df_results$symbol),
                         selected = "Fgf21", server=TRUE,options=list(maxItems=1))
    
    
    #update inputs for the heatplot
    updateNumericInput(session, inputId = "pValTile")
    updateSliderInput(session, inputId = "logFCTile")
    updateSliderInput(session, inputId = "nbTile")
    updateSliderInput(session, inputId = "nbGenes")
    
    #Create Forrest Plot
    output$forrestPlot <- renderPlot({
        
        #Select GEOSET's according to above input
        sex<-GEOSET_descr$GEOSET[GEOSET_descr$sex%in%input$GEOSET_sex]
        age<-GEOSET_descr$GEOSET[GEOSET_descr$age%in%input$GEOSET_age]
        diet<-GEOSET_descr$GEOSET[GEOSET_descr$diet_offspring%in%input$GEOSET_diet]
        
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
        df_select_gene<-volc_df[volc_df$symbol %in% input$gene_selectionV,]
        
        ggplot(volc_df, aes(x=logFC,y=-log10(P.Value), color=GEOSET))+
          geom_point(data=volc_df, aes(x=logFC,y=-log10(P.Value), color=GEOSET),
                     alpha=0.3, size=1.5)+
          geom_point(data=df_select_gene, aes(x=logFC,y=-log10(P.Value)),
                     color="black", shape=18, size=4)+
          theme_Publication()
    })
    
    #Expression rank plot
    #Define Unige color
    output$exprPlot <- renderPlot({
        color_unige<-c("#CF0063")
        
        #Import expression rank function
        R <- exp_rank()
        
        ggplot(reshape2::melt(R[input$gene_selection,,drop=FALSE])) + 
            geom_point(aes(x=Var2,y=value), shape=95, size=16,
                       color=color_unige) + 
            geom_col(aes(x=Var2, y=1), color="black", alpha=0.1)+
            theme_Publication() + 
            ylab("Gene Expression percentile") + 
            xlab("dataset") + 
            scale_y_continuous(labels=scales::percent,limits=c(0,1)) 
            #ggtitle(input$gene_selection)
        
    })
    
    #Create tile plot for display on the third row
    output$tilePlot<-renderPlot({
        ##Produce top gene list for tile plot
        #Select genes with specified criteria
        pre_tile <- df_results$symbol[abs(df_results$logFC)>input$logFCTile&
                                        df_results$P.Value<input$pValTile]
        genes_tiles <- names(which(table(pre_tile)>input$nbTile))
        
        #Producedf with top 40 genes
        df_genes_tile<-df_results[df_results$symbol%in%genes_tiles,]
        df_genes_tile<-df_genes_tile[order(abs(df_genes_tile$logFC), decreasing=TRUE),]
        genes_tile40<-unique(head(df_genes_tile, input$nbGenes)[,2])
        df_tile40<-df_results[df_results$symbol%in%genes_tile40,]
        
        #convert high logFC (smaller or larger than 2) to 2 and select appropriate columns
        df_tile40$logFC<-ifelse(df_tile40$logFC>2, 2,df_tile40$logFC)
        df_tile40$logFC<-ifelse(-2>df_tile40$logFC, -2,df_tile40$logFC)
        df_tile40<-df_tile40[,c(1:3)]
        
        #Produce the plot to display
        ggplot(df_tile40, aes(x = GEOSET, y = reorder(symbol, logFC), fill=logFC)) +   #reorder variables according to values
          scale_fill_gradient2(low="#2166ac", mid="#f7f7f7", high="#b2182b", limits=c(-2,2))+
            labs(x="GEO Expression Set", y="Gene Symbol", fill="logFC")+
            geom_tile(colour="white",size=0.1)+
            theme_Publication()
    })
    
    #Datatable with selected genes
    output$table <- renderDataTable({
        data <- df_results[df_results$symbol %in% input$gene_selection,
                           c(1:3, 6:8)]
        data$P.Value<-formatC(data$P.Value, format="E", digits=2)
        data$adj.P.Val<-formatC(data$adj.P.Val, format="E", digits=2)
        data
    })
    
    #Datatable with GEOSET description
    output$GEOSETtable <- renderDataTable({
      GEOSET_descr_tab
    })
    
}

shinyApp(ui = ui, server = server)