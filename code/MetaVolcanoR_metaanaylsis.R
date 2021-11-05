####################
#MetaVolcano Analsis
####################
#Metaanalysis of the dataset according to the MetaVolcanoR package
#Intitial analysis performed in June 2021 by Beat Moeckli

library(MetaVolcanoR)

#Import dataframe
df_results<-read.csv("df_results.csv")
list_results<-split(df_results, f=df_results$GEOSET)

#browseVignettes("MetaVolcanoR")

head(list_results[[9]])
length(list_results)

#Random effect model analysis
#The summary p-value stands for for the probability that the summary fold change
#is not different than Zero the perturbation ranking follows the topconfect approach
meta_degs_rem <- rem_mv(diffexp=list_results,
                        pcriteria="P.Value",
                        foldchangecol='logFC', 
                        genenamecol='symbol',
                        geneidcol=NULL,
                        collaps=TRUE,
                        llcol='CI.L',
                        rlcol='CI.R',
                        vcol=NULL, 
                        cvar=TRUE,
                        metathr=0.01,
                        jobname="MetaVolcano",
                        outputfolder=".", 
                        draw='HTML',
                        ncores=1)

#Gives the top perturbated genes, following the topconfect approach
head(meta_degs_rem@metaresult, 100)
#Draws MetaVolcano Plot
meta_degs_rem@MetaVolcano
#Draws forest plot for individual genes with with all the studies
draw_forest(remres=meta_degs_rem,
            gene="Fgf21",
            genecol="symbol", 
            foldchangecol="logFC",
            llcol="CI.L", 
            rlcol="CI.R",
            jobname="MetaVolcano",
            outputfolder=".",
            draw="HTML")

#assesses in how many of the datasets genes are differentially expressed.
#Criteria for the genes to be differentially expressed can be manually set.
meta_degs_vote <- votecount_mv(diffexp=list_results,
                               pcriteria="P.Value",
                               foldchangecol='logFC', 
                               genenamecol='symbol',
                               collaps=T,
                               pvalue=0.05,
                               foldchange=0, 
                               metathr=0.01)
meta_degs_vote@degfreq
head(meta_degs_vote@metaresult, 50000)
#Vulcano Graph showing consistency of differential expression
meta_degs_vote@MetaVolcano

#Analysis combining fold change and p-value combination. 
#"idx" is calculated as metafc*-log10(metap), this according to the topconfect method
meta_degs_comb <- combining_mv(diffexp=list_results,
                               pcriteria='P.Value', 
                               foldchangecol='logFC',
                               genenamecol='symbol',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.01, 
                               collaps=TRUE,
                               jobname="MetaVolcano",
                               outputfolder=".",
                               draw='HTML')

results_metavcomb<-head(meta_degs_comb@metaresult, 2500)
meta_degs_comb@MetaVolcano