#VISUALIZATION OF RESULTS FROM ALL ANALYSED DATASETS
#JUNE/JULY 2021, BEAT MOECKLI

#loading necessary packages
lapply(c("tidyverse", "RColorBrewer", "reshape", "data.table", "colorspace", "readxl"), require, character.only = TRUE)

mycolors<-c(brewer.pal(11, "RdBu"))
mycolors<-rev(mycolors)

#Import file containing all values for all analysed datasets
df_results<-read.csv("df_results.csv")

#Adult offspring datasets
adult_offspring <- c("GSE40903_CCvsOC",  "GSE40903_COvsOO",  "GSE123009_12w",  "GSE123009_28w",   
                     "GSE134976_CCvsOC", "GSE134976_COvsOO", "GSE44901_CCvsOC", "GSE44901_COvsOO",
                     "GSE46359_F", "GSE46359_M")
adult_offspring_CCvsOC<-c("GSE40903_CCvsOC",  "GSE123009_12w",  "GSE123009_28w",   
                          "GSE134976_CCvsOC", "GSE44901_CCvsOC")
adult_offspring_COvsOO<-c("GSE40903_COvsOO", "GSE134976_COvsOO", "GSE44901_COvsOO")
young_offspring<-c("GSE133767", "GSE46359_F", "GSE46359_M", "GSE62715")


###########################-
#Plot structure & functions
###########################


prep_melt<-function(gene_list){
  #Select only genes of interest from the whole dataset
  results_goi <-  dplyr::filter(df_results, symbol %in% gene_list)
  
  #Prepare table for further analysis
  MatrixData<-dplyr::select(results_goi, c(1,3,4))
  MatrixData<-as.data.table(MatrixData)
  MatrixData<-melt(MatrixData)
  MatrixData
}

tile_ggplot<-list(geom_tile(colour="white",size=0.1),
                  labs(x="GEO Expression Set", y="Gene Symbol", fill="logFC"))

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
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

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
  
}

###############################################
#Create graph with predefined genes of interest
###############################################

#Import a list of defined Genes of Interest
vec_goi<-read_excel("overview.xlsx", 2)%>%
  pull(Symbol) #Extract gene names in the symbol column of the overview file

#prepare Data for graph with above function
res_goi<-prep_melt(vec_goi)%>%dplyr::filter(GEOSET%in%adult_offspring_CCvsOC)

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_goi, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +  #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-3.1,3.1), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Genes of interest Expression Matrix")+tile_ggplot+
  theme_Publication()

#####################################################
#Create Matrix with top 5 expressed genes per GEO Set
#####################################################

#Filter top regulated genes
vec_topreg<-df_results%>%dplyr::filter(GEOSET%in%adult_offspring_CCvsOC)%>% #Filter out dataset from embryo data
  dplyr::filter(P.Value<=0.05)%>% #Filter out genes with a p-value above 0.05
  group_by(GEOSET)%>%top_n(5, abs(logFC))%>%  #select the top 5 genes per dataset
  pull(symbol)  #Pull symbol

#prepare Data for graph with above function
res_topreg<-prep_melt(vec_topreg)%>%dplyr::filter(GEOSET%in%adult_offspring_CCvsOC)

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_topreg, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +   #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-3.2,5.6), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Top 5 genes per dataset")+
  tile_ggplot+theme_Publication()


####################################################################
#Select genes with p<0.005, logFC>abs(0.2) and regulated >2 <GEOSETs
####################################################################

#Filter top regulated genes
vec_mul<-
  df_results%>%dplyr::filter(GEOSET!=c("GSE133767"))%>% #Filter out dataset from embryo data
  dplyr::filter(P.Value<=0.005)%>% #Filter out genes with a p-value above 0.05
  dplyr::filter(abs(logFC)>=0.2)%>%
  mutate(dir_up=case_when(logFC>0~1, logFC<0~0, TRUE~0))%>%
  mutate(dir_dn=case_when(logFC>0~0, logFC<0~1, TRUE~0))%>%
  group_by(symbol)%>%summarise(cons_up=sum(dir_up), cons_dn=sum(dir_dn), sum_logFC=sum(logFC))%>%
  dplyr::filter(cons_up>=3|cons_dn>=3)%>% dplyr::filter(abs(sum_logFC)>=1)%>%
  pull(symbol)

#prepare Data for graph with above function
res_mul<-prep_melt(vec_mul)%>%dplyr::filter(GEOSET!="GSE133767")

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_mul, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +   #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-2,2), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Genes of interest Expression Matrix")+
  tile_ggplot+ theme_Publication()

#####################################################################
#Select genes that are regulated in the same direction in (all) GEOSETS
#####################################################################

#Filter top regulated genes
vec_cons_genes<-
  df_results%>%dplyr::filter(GEOSET!=c("GSE133767"))%>% #Filter out dataset from embryo data
  dplyr::filter(P.Value<=0.5)%>% #Filter out genes with a p-value above 0.05
  dplyr::filter(abs(logFC)>=0.05)%>%
  mutate(dir=case_when(logFC>0~1, logFC<0~-1, TRUE~-1))%>%
  group_by(symbol)%>%summarise(cons=sum(dir))%>%
  dplyr::filter(abs(cons)>7)%>%
  pull(symbol)

#Filter top regulated genes
vec_cons_genes<-
  df_results%>%dplyr::filter(GEOSET!=c("GSE133767"))%>% #Filter out dataset from embryo data
  dplyr::filter(P.Value<=0.5)%>% #Filter out genes with a p-value above 0.05
  dplyr::filter(abs(logFC)>=0.05)%>%
  mutate(dir_up=case_when(logFC>0~1, logFC<0~0, TRUE~0))%>%
  mutate(dir_dn=case_when(logFC>0~0, logFC<0~1, TRUE~0))%>%
  group_by(symbol)%>%summarise(cons_up=sum(dir_up), cons_dn=sum(dir_dn), sum_logFC=sum(logFC))%>%
  dplyr::filter(cons_up>=8|cons_dn>=8)%>% dplyr::filter(abs(sum_logFC)>=1)%>%
  pull(symbol)

#prepare Data for graph
res_cons<-prep_melt(vec_cons_genes)%>%dplyr::filter(GEOSET!=c("GSE133767"))

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_cons, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +   #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-2,2), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Genes of interest Expression Matrix")+
  tile_ggplot+theme_Publication()

#####################################################################
#Select genes that are regulated in the same direction in >=4 GEOSETS
#plus selected on a sum of the logFC combined between all samples
#####################################################################

#Filter top regulated genes, tailored to include interesting candidates such as Fgf21
vec_cons_genes2<-
  df_results%>%dplyr::filter(GEOSET%in%adult_offspring)%>% #Filter out dataset from embryo data
  dplyr::filter(adj.P.Val<=0.9)%>% #Filter out genes with a p-value above 0.05
  dplyr::filter(abs(logFC)>=0.1)%>%
  mutate(dir_up=case_when(logFC>0~1, logFC<0~0, TRUE~0))%>%
  mutate(dir_dn=case_when(logFC>0~0, logFC<0~1, TRUE~0))%>%
  group_by(symbol)%>%summarise(cons_up=sum(dir_up), cons_dn=sum(dir_dn), sum_logFC=sum(logFC))%>%
  dplyr::filter(cons_up>=5|cons_dn>=5)%>% dplyr::filter(abs(sum_logFC)>=2)%>%
  pull(symbol)

#With Various filter options regarding conditions, CO vs OO
vec_cons_genes_COvsOO<-
  df_results%>%dplyr::filter(GEOSET%in%adult_offspring_COvsOO)%>% #Filter out datasets
  dplyr::filter(P.Value<=0.4)%>% #Filter out genes with a p-value above 0.05
  dplyr::filter(abs(logFC)>=0.15)%>%
  mutate(dir=case_when(logFC>0~1, logFC<0~-1, TRUE~0))%>%
  group_by(symbol)%>%summarise(cons=sum(dir), sum_logFC=sum(logFC))%>%
  dplyr::filter(abs(cons)>=2)%>% dplyr::filter(abs(sum_logFC)>=0.91)%>%
  pull(symbol)

#prepare Data for graph
res_cons2<-prep_melt(vec_cons_genes2)%>%
  dplyr::filter(GEOSET%in%adult_offspring)

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_cons2, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +   #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-3.1,3.1), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Genes of interest Expression Matrix")+
  tile_ggplot+theme_Publication()