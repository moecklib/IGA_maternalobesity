#ANALYSIS OF QPCR RESULTS FROM NOVEMBER 2021 (16 GENES IDENTIFIED DURING IGA)
#NOVEMBER 2021, BEAT MOECKLI

#loading necessary packages
lapply(c("tidyverse", "RColorBrewer", "reshape", "data.table", "colorspace", "readxl", 
         "ggpubr"), require, character.only = TRUE)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Raw results import and annotation####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Data file Impor, qPCR already analysed by the macro Excel document of the platform
qPCR_results_pre<-read.csv("data/qPCR_analysis/20211119_qPCR_RawResBefIntNorm.csv")

#Import of annotatation file 1 (cDNA list from Stephanie Lacotte),
#Cage number conlicts with the sampleKey2 file resolved manually
sampleKey1<-read.csv("data/qPCR_analysis/20211119_qPCR_sampleKey1.csv")%>%
  mutate(cage_ID=str_extract(Cage,"\\-.*"))%>%
  mutate(cage_ID=str_extract(cage_ID,"[:digit:]"))%>%
  transmute(cage=str_replace(Cage,"\\-.*",""), cage_ID=cage_ID, sample=sample)

#Import of annotatation file 1 (copy of no_DEN file with columns removed)
sampleKey2<-read.csv("data/qPCR_analysis/20211119_qPCR_sampleKey2.csv")%>%
  transmute(cage_ID=as.character(cage_ID),cage=str_replace(cage,"\\-",""),
            sex=sex, group=group, group2=group2)

#Join the different files together to obtain the necessary groups for grapsh(HFD vs ND, M_HFD vs M_ND)
qPCR_results<-qPCR_results_pre%>%left_join(sampleKey1)%>%
  left_join(sampleKey2, by = c("cage", "cage_ID"))%>%
  mutate(sample=as.character(sample))%>%
  relocate(where(is.numeric), .after = where(is.character))

#Create annotated outputfile, include Fgf21 manually and reread the file
#write.csv(qPCR_results, file="data/qPCR_analysis/20211119_qPCR_ResAnnotated.csv")
#qPCR_results<-read.csv(file="data/qPCR_analysis/20211119_qPCR_ResAnnotated.csv")

#Pivot longer to be able to display results in a facet wrap
long_qPCR_results<-qPCR_results%>%mutate(sample=as.character(sample), cage_ID=as.character(cage_ID))%>%
  pivot_longer(where(is.numeric), names_to = "gene", values_to = "Ct")%>%
  group_by(gene, group)%>%
  dplyr::mutate(mean=mean(Ct))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Plot structure & functions####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

theme_Publication <- function(base_size=32, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(t=0,r=0,b=25,l=0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x  = element_text(), 
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


qPCR_style<-list(geom_boxplot(size=1.5),
                    geom_beeswarm(color="black", size=5, alpha=0.8, shape=17),
                    geom_signif(comparisons = list(c("F_HFD", "F_ND")), test = wilcox.test, textsize = 11,
                                size=1),
                    geom_signif(comparisons = list(c("M_HFD", "M_ND")), test = wilcox.test, textsize = 11,
                                size=1),
                    theme_Publication(),
                    theme(legend.position="none")
)

mycolors2<-c(brewer.pal(8, "Paired")[c(8,7,2,1)])

#Function to create qPCR plots
qPCR_plot<-function(gene){
  ggplot(data=long_qPCR_results[long_qPCR_results$gene%in%gene,], 
         aes(x=group, y=Ct))+
    ggtitle(gene)+
    labs(y=paste("Relative",gene, "Expression"), x=NULL)+
    qPCR_style
}

#Save the plot in a standardized format
save_plot<-function(plot){
  ggsave(
    
    filename=paste(plot,".tiff"),
    device="tiff",
    width=12,
    height=10,
    # units="mm",
    path="output/graphs/Figure_5",
    dpi = "retina")
}


#*#*#*#*#*#*#*#*#*#*#*#
#Plots for figure 5####
#*#*#*#*#*#*#*#*#*#*#*#

#Expression data in facet plot to visualize all experiments together
ggplot(data=long_qPCR_results, aes(x=group, y=Ct, fill=group))+geom_boxplot(show.legend = FALSE)+
  facet_wrap(~gene, scales="free")+
  stat_compare_means(comparisons=list(c("M_HFD", "M_ND")),label="p.format", hide.ns = TRUE)+
  stat_compare_means(comparisons=list(c("F_HFD", "F_ND")),label="p.format", hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                     trans="log2")+
  scale_fill_manual(values=c(mycolors2))+
  theme_Publication()+labs(y= NULL, x=NULL)+ ggtitle("Expression data in No DEN animals")

#Create plots for the figure
qPCR_plot("Ppara")
qPCR_plot("Fgf21")


#Execute the function to save the plot
save_plot("Ppara")

