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

#Import main file with all the results & description GEOSETS
df_results<-read.csv("output/df_results.csv")
GEO_descr<-read.csv("data/description_GEOSET.csv")

##Aldh2 analysis
#Filter all studies with the expression of Aldh2
Aldh2_pre<-df_results[df_results$symbol=="Aldh2",]
Aldh2<-merge(x=Aldh2_pre, y=GEO_descr[,c("GEOSET", "Author", "Year", "Offspring_diet",
                                         "Sex", "Age")],
             by="GEOSET", all.x=TRUE)%>%
  arrange(Year)

#Generation of meta object for further analysis
Aldh2.meta<-metagen(TE=logFC,
        studlab = Author,
        data=Aldh2,
        lower=CI.L,
        upper=CI.R,
        pval=adj.P.Val)

colnames(Aldh2)
?metagen

#Creation of plot
forest.meta(Aldh2.meta,
            layout="JAMA")
forest.meta(Aldh2.meta,
            fixed=F,
            print.tau2=F,
            col.diamond = "tomato1",
            leftcols =c("studlab", "Year", "Offspring_diet",
                      "Sex", "Age"))

?forest.meta

##Acsl1 analysis
#Filter all studies with the expression of Acsl1
Acsl1_pre<-df_results[df_results$symbol=="Acsl1",]
Acsl1<-merge(x=Acsl1_pre, y=GEO_descr[,c("GEOSET", "Author", "Year", "Offspring_diet",
                                         "Sex", "Age")],
             by="GEOSET", all.x=TRUE)%>%
  arrange(Year)

#Generation of meta object for further analysis
Acsl1.meta<-metagen(TE=logFC,
                    studlab = Author,
                    data=Acsl1,
                    lower=CI.L,
                    upper=CI.R,
                    pval=adj.P.Val)

forest.meta(Acsl1.meta,
            layout="JAMA")
forest.meta(Acsl1.meta,
            fixed=F,
            print.tau2=F,
            col.diamond = "tomato1",
            leftcols =c("studlab", "Year", "Offspring_diet",
                        "Sex", "Age"))

comb.meta<-metamerge(Acsl1.meta, Aldh2.meta,
                     text.pooled1 = "Acsl1",
                     text.pooled2 = "Aldh2")
forest.meta(comb.meta,
            print.tau2=F,
            col.diamond = "tomato1",
            leftcols =c("studlab", "Year", "Offspring_diet",
                        "Sex", "Age"))


            fixed=F,
            print.tau2=F,
            col.diamond = "tomato1",
            leftcols =c("studlab", "Year", "Offspring_diet",
                        "Sex", "Age"))
?metabind
?metamerge
