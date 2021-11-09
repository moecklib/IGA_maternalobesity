#GENERATION OF TOP DIFFERENTIALLY EXPRESSED GENES
#NOVEMBER 2021, BEAT MOECKLI, JULIEN PRADOS

#loading necessary packages
lapply(c("tidyverse", "RColorBrewer", "reshape", "data.table", "colorspace", "readxl", "EnsDb.Mmusculus.v79",
         "clusterProfiler"), require, character.only = TRUE)

mycolors<-c(brewer.pal(11, "RdBu"))
mycolors<-rev(mycolors)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Plot structure & functions####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

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

#*#*#*#*#*#*#*#*#*#*#*#*
#Import base dataset####
#*#*#*#*#*#*#*#*#*#*#*#*
#Import file containing all values for all analysed datasets
df_results<-read.csv("df_results.csv")

#Remove datasets with fetal or blastcyst samples
df_results<-read.csv("df_results.csv")%>%
  dplyr::filter(!GEOSET%in%c("GSE133767", "GSE62715")) #Filter datasets from fetal and blastocyst offspring

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Perform the selection of genes according to the below steps####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
top_DEG<-local({

#1. Define selection criteria for the datasets
thrLFC<-  0.1
thrnegLFC<--0.1
thrpVal<- 0.05

#2. Add additional variables to the datasets for later selection
top_DEGpreP<-df_results%>%group_by(GEOSET)%>%
  mutate(select=sum(abs(logFC)>thrLFC & P.Value<thrpVal))%>%
  add_count(name="total")%>%
  mutate(FDR=select/(select+total))%>%
  mutate(dir_up=case_when(logFC>thrLFC & P.Value<thrpVal~1, TRUE~0))%>%
  mutate(dir_dn=case_when(logFC<thrnegLFC & P.Value<thrpVal~1, TRUE~0))%>%
  mutate(DEG=case_when(abs(logFC)>thrLFC & P.Value<thrpVal~FDR, TRUE~1-FDR))%>%
  ungroup()%>%
  group_by(symbol)%>%
  mutate(cons_up=sum(dir_up),cons_dn=sum(dir_dn), 
         mean_logFC=mean(abs(logFC)), FDR_comb=prod(DEG), 
         Fisher_P.Val=pchisq((sum(log(P.Value))*-2), df=length(P.Value)*2, lower.tail=F))%>%
  ungroup()

#3. Perform the selection of the genes with the following criteria
#logFC>0.1, P.Val<0.05, sumLogFC>1, dysregulated in more than 3 conditions
top_DEG<-top_DEGpreP%>%
  #dplyr::filter(abs(logFC)>thrLFC & P.Value<thrpVal)%>%
  group_by(symbol)%>%
  summarise(cons_up=unique(cons_up),cons_dn=unique(cons_dn), 
            mean_logFC=unique(mean_logFC), FDR_comb=unique(FDR_comb), 
            Fisher_P.Val=unique(Fisher_P.Val))%>%
  dplyr::filter(cons_up>=3 |cons_dn>=3)%>%
  dplyr::filter(mean_logFC>=0.1)%>%
  arrange(desc(cons_up, abs(cons_dn), FDR_comb, mean_logFC))
})

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Add pathway information to genelist####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

top_DEG_annocleaned<-local({
#Get the wikipathways annotation file
wp2gene<-read.gmt.wp("data/wikipathways-20210610-gmt-Mus_musculus.gmt")

#obtain a list with gene names and EntrezID, command "columns(EnsDb.Mmusculus.v79)" shows available colnames
#Possibility to add entrez ID to gene names
gene_names<-AnnotationDbi::select(EnsDb.Mmusculus.v79, keys(EnsDb.Mmusculus.v79), 
                                  columns=c("ENTREZID", "GENENAME"))

top_DEG$entrezID<-gene_names$ENTREZID[match(top_DEG$symbol, gene_names$GENENAME)] #Add additional variable with EntrezID
top_DEG$entrezID<-as.character(top_DEG$entrezID)    #change EntrezID to character variable for left_join

#Left join the pathway information from wikipathways with the results table
top_DEG_anno<-top_DEG%>%left_join(wp2gene%>%dplyr::select(name, gene),
                                  by=c("entrezID"="gene"))%>%
  group_by(symbol,cons_up, cons_dn, mean_logFC, FDR_comb, Fisher_P.Val)%>%
  summarize(pathway=paste(name, collapse = ", "))

test<-top_DEG%>%left_join(wp2gene%>%dplyr::select(name, gene),
                                  by=c("entrezID"="gene"))%>%
  group_by(name)%>%
  summarize(pathway=n())%>%
  arrange(desc(pathway))

#Filter all rows with no pathway attributed
top_DEG_annocleaned<-top_DEG_anno%>%dplyr::filter(pathway!="NA")%>%
  mutate(sum_cons=cons_up+cons_dn, .after=cons_dn)%>%
  arrange(desc(sum_cons), FDR_comb, Fisher_P.Val)
})

#Write the results file to the output folder
#write.csv(top_DEG_annocleaned, file="output/TopDE_Genes.csv")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Visualisations of top dereg genes####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Filter top regulated genes identified according to the below criteria 
#(logFC>0.1, P.Val<0.05, sumLogFC>1, more than 3 conditions)
vec_top_DEG<-top_DEG%>%pull(symbol)

#prepare Data for graph
res_top_DEG<-prep_melt(vec_top_DEG)

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_top_DEG, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +   #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-2,2), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Genes of interest Expression Matrix")+
  tile_ggplot+theme_Publication()

#Only plot genes with relevant pathway information 
#Genes contained in the top_DEG_anno_cleaned dataframe
vec_top_DEG_path<-top_DEG_anno_cleaned%>%
  pull(symbol)

#prepare Data for graph
res_top_DEG_path<-prep_melt(vec_top_DEG_path)

#Plot data with geom_tile function, facet grid according to species)
ggplot(res_top_DEG_path, aes(x = GEOSET, y = reorder(symbol, value), fill=value)) +   #reorder variables according to values
  scale_fill_continuous_divergingx(palette="RdBu", limits=c(-2,2), rev=TRUE, mid = 0, l3 = 0, p1 = .2, p2 = .6, p3=0.6, p4=0.8) +
  labs(title="Genes of interest Expression Matrix")+
  tile_ggplot+theme_Publication()