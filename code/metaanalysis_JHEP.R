#GENERATION OF META ANALYSIS FORREST PLOT FOR Acsl1 and Aldh2
#January 2022, BEAT MOECKLI
#Analysis for lette to the editors in JHEP in response to Sun et al, J Hepatol., 2020


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Loading necessary packages####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#loading necessary packages
lapply(c("tidyverse", "meta"), require, character.only = TRUE)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Analysis with the meta package and generation of plots####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#Import main file with all the results
df_results<-read.csv("output/df_results.csv")

##Aldh2 analysis
#Filter all studies with the expression of Aldh2
Aldh2<-df_results[df_results$symbol=="Aldh2",]

#Generation of meta object for further analysis
Aldh2.meta<-metagen(TE=logFC,
        studlab = GEOSET,
        data=Aldh2,
        lower=CI.L,
        upper=CI.R)

#Creation of plot
forest.meta(Aldh2.meta,
            layout="JAMA")

##Acsl1 analysis
#Filter all studies with the expression of Acsl1
Acsl1<-df_results[df_results$symbol=="Acsl1",]

Acsl1.meta<-metagen(TE=logFC,
                    studlab = GEOSET,
                    data=Acsl1,
                    lower=CI.L,
                    upper=CI.R)

forest.meta(Acsl1.meta,
            layout="JAMA")