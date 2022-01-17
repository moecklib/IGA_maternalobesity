#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Plot structure & functions####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#loading necessary packages
lapply(c("tidyverse", "RColorBrewer", "reshape", "data.table", "colorspace",
         "ggbeeswarm", "ggpubr"), require, character.only = TRUE)


#Load colors
mycolors_groups<-c(brewer.pal(8, "Paired")[c(8,7,2,1)])

#Load publication theme
theme_Publication <- function(base_size=32, base_family="sans") {
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
            axis.text.x  = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.box = "vertical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

boxplot_style<-list(geom_boxplot(size=1.5),
                        geom_beeswarm(color="black", size=5, alpha=0.8, shape=17),
                        geom_signif(comparisons = list(c("F_HFD", "F_ND")), test = wilcox.test, textsize = 11,
                                    size=1),
                        geom_signif(comparisons = list(c("M_HFD", "M_ND")), test = wilcox.test, textsize = 11,
                                    size=1),
                        scale_fill_manual(values=mycolors_groups),
                        theme_Publication(),
                    theme(legend.position="none")
)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Import datafiles####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

histo_Tihy<-read.csv("data/No_DEN_Tihy2.csv")
no_DEN<-read.csv("data/No_DEN_1.csv")
weight_offspring_pre<-read.csv("data/weight_offspring_NoDEN.csv")

#Create long format file for the weight plot
weight_offspring<-weight_offspring_pre%>%
  pivot_longer(where(is.numeric), names_to = "week", values_to = "weight")%>%
  mutate(week=as.numeric(gsub(".*?([0-9]+).*", "\\1", week)))

#*#*#*#*#*#*#*#*#
#Create plots####
#*#*#*#*#*#*#*#*#

#Boxplot graph with all 40 timepoints four groups per timepoint. 
#Not used for the final figure because to busy and lacking M-ND animals for the later timepoints
ggplot(data=weight_offspring, aes(x=group, y=weight, fill=group))+
  geom_boxplot()+
  geom_beeswarm(color="black", size=1, alpha=0.6)+
  facet_wrap(~week, nrow=1, switch = "x")+
  stat_compare_means(comparisons=list(c("M_HFD", "M_ND")),label="p.signif", hide.ns = TRUE)+
  stat_compare_means(comparisons=list(c("F_HFD", "F_ND")),label="p.signif", hide.ns = TRUE)+
  scale_fill_manual(values=mycolors_groups)+
  labs(y= "Weight", x="Week after weaning")+
  theme_Publication()+
  ggtitle("Weight curve cohousing offspring")

#Weight curve graph for the first 21 weeks.
#Not used since esthetically less pleasing and lacking M-ND animals for the later timepoints
weight_offspring%>%dplyr::filter(week<21)%>%
ggplot(aes(x=week, y=weight, fill=group, color=group))+
  #geom_beeswarm(size=1, alpha=0.6)+
  geom_smooth()+
  scale_fill_manual(values=mycolors_groups)+
  scale_color_manual(values=mycolors_groups)+
  labs(y= "Weight [g]", x="Week after weaning")+
  theme_Publication()

#Weight plot final with legend and x-axis labels
weight_offspring%>%dplyr::filter(week%in%c(4,18))%>%
  ggplot(aes(x=group, y=weight, fill=group))+
  geom_boxplot(size=1.5)+
  geom_beeswarm(color="black", size=5, alpha=0.8, shape=17)+
  facet_wrap(~week, nrow=1, switch = "x")+
  geom_signif(comparisons = list(c("F_HFD", "F_ND")), test = wilcox.test, textsize = 9,
              size=1)+
  geom_signif(comparisons = list(c("M_HFD", "M_ND")), test = wilcox.test, textsize = 9,
              size=1)+
  scale_fill_manual(values=mycolors_groups)+
  labs(y= "Weight [g]", x="Week of life")+
  theme_Publication()+
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 0.8), 
        axis.ticks.x=element_blank())+
  guides(shape=FALSE)+
  ggtitle("Weight 4 & 18 weeks")+
  scale_y_continuous(limits=c(12, 34), expand = c(0.1, 0))    #Make the p-value labels visible

#ALT plot
ggplot(data= no_DEN, aes(x = group, y = ALT_24w, fill=group))+
  boxplot_style+
  scale_y_log10(limits = c(0.8e1, 1.5e2), expand = c(0, 0))+
  labs(y= "log10( ALT conc [IU/l] )", x=NULL)+
  ggtitle("Alanine transaminase")


#OGTT plot
ggplot(data= no_DEN, aes(x = group, y = OGTT_24w, fill=group))+
  boxplot_style+
  labs(y= "AUC [mmol/min]", x=NULL)+
  ggtitle("Oral Glucose Tolerance")


#NAFLD Barplots
ggplot(data= no_DEN, aes(x = group, fill=Bedossa_et_al.))+
  geom_bar(position="fill", color="white")+
  geom_text(data = no_DEN %>% group_by(group, Bedossa_et_al.) %>% tally() %>% #Add labels to barplot
              mutate(p = n / sum(n)) %>% ungroup(),
            aes(y = p, label = scales::percent(p)),
            position = position_stack(vjust = 0.5),
            show.legend = FALSE,
            size=8)+
  labs(y= "Proportion Histology", x=NULL)+
  scale_fill_manual(values = c("No NAFLD" = "#cccccc", "NAFLD" = "#fd8f24", "NASH" = "#c03728"))+
  ggtitle("Proportion of NAFLD or NASH")+
  theme_Publication()+
  theme(legend.title = element_blank())


#Fibrosis Boxplot
ggplot(data= no_DEN, aes(x = group, y = rel_fibrosis_Qpath, fill=group))+
  boxplot_style+
  scale_y_continuous(trans='log2',
                     expand = expansion(mult = c(0, 0.15)))+
  labs(y= "Surface fibrosis [%]", x=NULL)+
  ggtitle("Fibrosis")

#Steatosis Boxplot
ggplot(data= no_DEN, aes(x = group, y = MT2_ratio_surface, fill=group))+
  boxplot_style+
  scale_y_continuous(trans='log2', #log2 transformation of y-scale
                     expand = expansion(mult = c(0, 0.15)))+ #Add additional space
  labs(y= "Surface steatosis [%]", x=NULL)+
  ggtitle("Steatosis")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Test for steatosis calculation with maximum value of the surface ratio####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#Import data file from Mathieu
histo_Tihy<-read.csv("data/No_DEN_Tihy2.csv")

#Extract higher Ratio Surface of the two measurements
histo_tihy_max<-histo_Tihy%>%group_by(cage, cage_ID, group, group2)%>%
  dplyr::summarise(max_surfRatio=max(RatioSurfaceSurSurfaceTotale),
                max_vacRatio=max(RatioCelluleSurVacuoles))

#Display the corresponding graph
ggplot(data= histo_tihy_max, aes(x = group2, y = max_surfRatio, fill=group2))+
  geom_boxplot(size=1.5)+
  geom_beeswarm(color="black", size=5, alpha=0.8, shape=17)+
  geom_signif(comparisons = list(c("HFD", "ND")), test = wilcox.test, textsize = 11,
            size=1)+
  scale_y_continuous(trans='log2')+
  theme_Publication()+
  theme(legend.position="none")+
  labs(y= "log2 (Relative surface of steatosis [%])", x=NULL)+
  ggtitle("Steatosis")
