#Standardized analysis of MA expression set for later analysis in MetaVolcano
#Project to analyze the differential expression between offspring of high fat diet mothers and
#offspring of normal diet mothers. Different species are included.
#In this part of the project the datafile containing the different GEO sets is "created"
#Code developed by Beat Moeckli, University of Geneva, June 2021

#General analysis parameters used for all datasets:
#-Cutoff for low expressed genes 20th percentile, at least two samples with expression

#BiocManager::install(EnsDb.Mmusculus.v79)
#install.packages("BiocManager")
#install.packages("limma")

#loading necessary packages
lapply(c("GEOquery", "forcats", "stringr", "BiocManager", "ggrepel", "survminer", "limma", "pheatmap", 
         "org.Hs.eg.db", "RColorBrewer", "dplyr", "tidyr", "tidyverse", "annotation", 
         "EnsDb.Mmusculus.v79"), require, character.only = TRUE)

# sessionInfo()


#Load some Functions####

#Function to produce fit. Takes GSE as input and requires design and contrasts
produce_fit<-function(gse_xy){
  gse<-gse_xy
  cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
  is_expressed <- exprs(gse) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
  keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
  # #table(keep)                             # check how many genes are removed / retained.
  gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
  fit <- lmFit(exprs(gse), design)
  fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
  fit2 <- eBayes(fit2)
  fit2
}

#Clean up full result sheet. Takes full results (produced by topTable) as input and returns cleaned dataframe.
cleanup_results<-function(full_results){
  full_results%>%
    na.omit(cols=c("symbol", "logFC"))%>%  #omit rows/probes with no symbol assigned
    group_by(symbol)%>%
    arrange(desc(AveExpr))%>%
    slice_head(n = 1)  #select top signal probe for each gene
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
##Load all individual datasets into functions####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#*#*#*#*#*
#GSE133767
#*#*#*#*#*

GSE133767<-function(){
  #Download expression dataset
  gse_133767 <- getGEO(filename = "data/GSE133767_series_matrix.txt.gz")
  
  #Obtain annotation for samples
  anno_133767<-local({
    #Determine annotation of samples with the data stored under feature Data
    anno <- fData(gse_133767)%>%                     #Obtain a list with information for each probe
      dplyr::select(ID,gene_assignment)       #Select columns that are of interest
    
    #Select ENS items from the column gene information
    anno <- str_split(anno$gene_assignment, boundary("word"))%>%          #Split up the column in different strings/words
      enframe()%>%                                                        #Recreate a dataframe with a column-list
      cbind(anno)%>%                                                      #Combine with the original df to not use probe ID
      unnest(value)%>%dplyr::select(value, ID)%>%                         #Create a long format table storing all elements from "gene_id"
      dplyr::filter(str_detect(value, "^ENS")) 
    
    #match ENSMUSTXID to the corresponding gene symbols, command "columns(EnsDb.Mmusculus.v79)" shows available colnames
    gene_names<-select(EnsDb.Mmusculus.v79, keys(EnsDb.Mmusculus.v79), c("TXID", "GENEID", "GENENAME"))
    
    anno$symbol<-gene_names$GENENAME[match(anno$value, gene_names$TXID)]
    anno<-dplyr::select(anno, ID, symbol)%>%distinct() #keeps only the first row if several row of a unique combination
    anno$ID<-as.character(anno$ID)
    anno
  })
  
  #Perform LIMMA for GSE133767 and produce fit
  fit_133767<-local({
    #Select the liver samples, the sample ID's have been identified with the GEO2R tool
    samples<-c((1:8),(17:24))   #Create a vector that corresponds to the liver samples
    gse <- gse_133767[ ,samples]               #selects only the liver samples
    
    # exprs get the expression levels as a data frame and get the distribution
    summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    #exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    ##Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             
    
    ## Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- dplyr::select(sampleInfo, characteristics_ch1)%>%
      dplyr::rename(F0_diet=characteristics_ch1)                                #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-transmute(sampleInfo, dietcombo=case_when(F0_diet=="diet of parents: Control Diet (CD) (lean parents)"~"CC",
                                                          F0_diet=="diet of parents: High fat diet (HFD) (obese mother)"~"OC",
                                                          TRUE~"MOC"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC","OC")         #Rename Col names to make it prettier
    
    
    ##Filter low expressed genes from analysis
    # cutoff <- quantile(exprs(gse), 0.2, na.rm=TRUE)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # table(keep)                             # check how many genes are removed / retained.
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    ##Determine the expression level in each pre-specified group
    fit <- lmFit(exprs(gse), design)
    
    #specify contrasts
    contrasts <- makeContrasts(CC - OC, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit2 <- eBayes(fit2)
    
    topTable<-topTable(fit2)              #Get first look at the results                                
    table(decideTests(fit2))              #To see how many genes are deferentially expressed
    fit2
  })
  
  #Create results and clean up result table
  results_133767<-local({
    #Write a table with the full results. Specify the number of rows (inf=infinite)
    full_results <- topTable(fit_133767, number=Inf, confint = T)%>%    #Include ConfInt
      tibble::rownames_to_column("ID")
    
    #Insert the new annotation into the full results table. Join the two tables by the variable ID
    full_results<-left_join(anno_133767, full_results, by=c("ID" = "ID"))
    
    #Clean up full result sheet
    full_results_clean<-cleanup_results(full_results)
    full_results_clean
  })
  
  #Create a list with the result files
  GSE133767<-list(GSE133767=results_133767)
}

#**********************************

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#GSE40903_CCvsOC&GSE40903_OOvsCO
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

GSE40903<-function(){
  
  #Download expression dataset
  gse_40903 <- getGEO(filename = "data/GSE40903_series_matrix.txt.gz")
  
  #Create fit including both conditions (F1HFD and F1ND)
  fit_40903<-local({
    
    #Select the liver samples, the sample ID's have been identified with the GEO2R tool
    samples<-seq(81, length.out = 36)   #Create a vector that corresponds to the liver samples
    gse <- gse_40903[ ,samples]               #selects only the liver samples
    
    #exprs get the expression levels as a data frame and get the distribution
    #summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    ##Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             
    
    ## Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- dplyr::select(sampleInfo, characteristics_ch1.3, characteristics_ch1.4)%>%
      dplyr::rename(F0_diet=characteristics_ch1.3, F1_diet=characteristics_ch1.4)              #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-transmute(sampleInfo, dietcombo=case_when(F0_diet=="maternal diet: LF"& F1_diet=="diet: LF"~"CC",
                                                          F0_diet=="maternal diet: HF"& F1_diet=="diet: LF"~"OC",
                                                          F0_diet=="maternal diet: LF"& F1_diet=="diet: HF"~"CO",
                                                          F0_diet=="maternal diet: HF"& F1_diet=="diet: HF"~"OO",
                                                          TRUE~"else"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC", "OC", "CO", "OO")         #Rename Col names to make it prettier
    
    ##Filter low expressed genes from analysis
    # cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # #table(keep)                             # check how many genes are removed / retained.
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    ##Determine the expression level in each pre-specified group
    fit <- lmFit(exprs(gse), design)
    
    #specify contrasts
    contrasts <- makeContrasts(CC - OC, CO-OO, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit2 <- eBayes(fit2)
    
    anno <- fData(gse)                    #Obtain a list with information for each probe
    anno <- dplyr::select(anno,Symbol)    #Select columns that are of interest
    
    fit2$genes <- anno                    #Insert this new information into the table  with the results
    fit2
  })
  
  #Extract results with topTable and clean up omitting NA only selecting one sample per gene
  results_40903_CCvsOC<- local({
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results_F1ND <- topTable(fit_40903, coef=1, number=Inf, confint=T)%>%   #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")                                          #ProbeID to a column rather than Row name
    
    #Clean up full result sheet
    full_results_clean_F1ND<-full_results_F1ND%>%
      dplyr::rename(symbol=Symbol)%>% #rename "Symbol" to match other datasets
      cleanup_results()
    full_results_clean_F1ND
  })
  
  results_40903_COvsOO<- local({
    full_results_F1HFD <- topTable(fit_40903, coef=2, number=Inf, confint=T)%>%    #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")
    
    #Clean up full result sheet
    full_results_clean_F1HFD<-full_results_F1HFD%>%
      dplyr::rename(symbol=Symbol)%>%  #rename "Symbol" to match other datasets
      cleanup_results()
    full_results_clean_F1HFD
    
  })
  
  #Create a list with the TWO new result files
  GSE40903<-list(GSE40903_CCvsOC=results_40903_CCvsOC, 
                  GSE40903_COvsOO=results_40903_COvsOO)
}

#**********************************

#*#*#*#*#*#*#*
#GSE123009_12w
#*#*#*#*#*#*#*

GSE123009_12w<-function(){
  #Download expression dataset
  gse_123009 <- getGEO(filename = "data/GSE123009_series_matrix.txt.gz")
  
  fit_123009_12w<-local({
    #Select the liver samples, the sample ID's have been identified with the GEO2R tool
    XIIw_samples<-c(c(1:7),c(17:23))        #Create a vector that corresponds to the liver and 12w samples
    gse <- gse_123009[ , XIIw_samples]       #selects only the liver samples
    
    #exprs get the expression levels as a data frame and get the distribution
    #summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    #Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             
    
    # Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- dplyr::select(sampleInfo, characteristics_ch1)%>%
      dplyr::rename(F0_diet=characteristics_ch1)              #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-transmute(sampleInfo, dietcombo=case_when(F0_diet=="prenatal diet: Low fat"~"CC",
                                                          F0_diet=="prenatal diet: High fat"~"OC",
                                                          TRUE~"NA"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC", "OC")         #Rename Col names to make it prettier
    
    ##Filter low expressed genes from analysis
    # cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    ##Determine the expression level in each pre-specified group
    fit <- lmFit(exprs(gse), design)
    
    #specify contrasts
    contrasts <- makeContrasts(CC - OC, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit2 <- eBayes(fit2)
    fit_123009_12w<-fit2
  })
  
  #Create output with annotation
  results_123009_12w<-local({
    
    #Determine annotation of samples with the data stored under 
    anno <- fData(gse_123009)                    #Obtain a list with information for each probe
    gene_names<-select(EnsDb.Mmusculus.v79, keys(EnsDb.Mmusculus.v79), 
                       columns=c("ENTREZID", "SYMBOL")) #list with information about used genes
    anno$symbol<-gene_names$SYMBOL[match(anno$ENTREZ_GENE_ID, gene_names$ENTREZID)] #match gene symbol and ENTREZ ID
    anno <- dplyr::select(anno,c(symbol))    #Select columns that are of interest (basically only symbol)
    
    #Insert this new information into the table  with the results
    fit_123009_12w$genes <- anno
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results<- topTable(fit_123009_12w, number=Inf, confint = T)%>%    #Specify which contrast you want to look if several
      tibble::rownames_to_column("ID")              #ProbeID to a column rather than Row name
    
    #Clean up full result sheet
    full_results_clean<-cleanup_results(full_results)
    
  })
  
  GSE123009_12w<-list(GSE123009_12w=results_123009_12w)
}

#**********************************

#*#*#*#*#*#*#*
#GSE123009_28w
#*#*#*#*#*#*#*

GSE123009_28w<-function(){
  #Download expression dataset
  gse_123009 <- getGEO(filename = "data/GSE123009_series_matrix.txt.gz")
  
  fit_123009_28w<-local({
    #Select the liver samples, the sample ID's have been identified with the GEO2R tool
    XXVIIIw_samples<-c(c(8:16),c(24:32))                     #Create a vector that corresponds to the liver samples
    gse <- gse_123009[ , XXVIIIw_samples]               #selects only the liver samples
    
    #exprs get the expression levels as a data frame and get the distribution
    summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    #Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             
    
    # Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- dplyr::select(sampleInfo, characteristics_ch1)%>%
      dplyr::rename(F0_diet=characteristics_ch1)              #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-transmute(sampleInfo, dietcombo=case_when(F0_diet=="prenatal diet: Low fat"~"CC",
                                                          F0_diet=="prenatal diet: High fat"~"OC",
                                                          TRUE~"NA"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC", "OC")         #Rename Col names to make it prettier
    
    ##Filter low expressed genes from analysis
    # cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    ##Determine the expression level in each pre-specified group
    fit <- lmFit(exprs(gse), design)
    
    #specify contrasts
    contrasts <- makeContrasts(CC - OC, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit2 <- eBayes(fit2)
    fit2
  })
  
  #Create output with annotation
  results_123009_28w<-local({
    
    #Determine annotation of samples with the data stored under 
    anno <- fData(gse_123009)                    #Obtain a list with information for each probe
    gene_names<-ensembldb::select(EnsDb.Mmusculus.v79, keys(EnsDb.Mmusculus.v79), 
                                  columns=c("ENTREZID", "SYMBOL")) #list with information about used genes
    anno$symbol<-gene_names$SYMBOL[match(anno$ENTREZ_GENE_ID, gene_names$ENTREZID)] #match gene symbol and ENTREZ ID
    anno <- dplyr::select(anno,c(symbol))    #Select columns that are of interest (basically only symbol)
    
    #Insert this new information into the table  with the results
    fit_123009_28w$genes <- anno
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results<- topTable(fit_123009_28w, number=Inf, confint = T)%>%    #Specify which contrast you want to look if several
      tibble::rownames_to_column("ID")              #ProbeID to a column rather than Row name
    
    #Clean up full result sheet
    full_results_clean<-full_results%>%
      cleanup_results()
    full_results_clean
    
  })
  
  GSE123009_28w<-list(GSE123009_28w=results_123009_28w)
}

#**********************************

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#GSE134976_CCvsOC&GSE134976_OOvsCO
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

GSE134976<-function(){
  
  #Import and transform expression data into an ExpressionSet
  gse_134976<-local({
    #Import file containing some sample and analysis information but no expression values
    gse_134976_pdata <- getGEO(filename = "data/GSE134976_series_matrix.txt.gz")
    #Adapt samples to correspond to conditions of interest
    samples<-c((1:5),(7:12))   #Create a vector that corresponds to the liver samples
    gse_134976_pdata <- gse_134976_pdata[ ,samples]
    
    #Import matrix containing expression data. Per above description the file contains " include RPKM values for each Sample"
    GSE134976_cpm <- read.table("data/GSE134976_cpm_all_groups_rerun_glm.txt.gz", header=TRUE, sep="\t")
    
    #Construct the Biobase ExpressionSet with the given elements
    gse_134976<-ExpressionSet(data.matrix(GSE134976_cpm)) #enter expression data in matrix form into the GSE
    pData(gse_134976)<-pData(gse_134976_pdata) #enter pData into the newly constructed expressionset
    
    #Enter fData with sample information as needed
    gene_names<-ensembldb::select(EnsDb.Mmusculus.v79, 
                                  keys(EnsDb.Mmusculus.v79), c("GENENAME", "GENEID")) #Obtain gene names for fData
    ?select
    rownames<-rownames(gse_134976)  #create vector with genenames, GENEID's
    dimnames<-c(list(rownames), c("GENEID")) #create dimension names of matrix (needs to be a list with two elements)
    fData<-matrix(rownames, ncol=1, dimnames=dimnames) #create matrix with two columns and GENEIDs to match
    fData<-fData%>%data.frame()%>% #convert matrix to a dataframe
      left_join(gene_names, by="GENEID")%>%  #left_join thecreated datamatrix with the above genenames (pipe, first element is fData)
      mutate(GENE_ID=GENEID)%>%   #Add a duplicate column for GENEID for safety
      column_to_rownames(var="GENEID")  #Change the rowname to correspond to the rowname of the exprs
    fData(gse_134976)<-fData            #Enter the fData into the eSet
    colnames(exprs(gse_134976))<-rownames(pData(gse_134976))
    
    #summary(exprs(gse_134976))    #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    exprs(gse_134976) <- log2(exprs(gse_134976))  #log2 transformation of the expression values
    #boxplot(exprs(gse_134976),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    gse_134976
  })
  
  #Produce fit
  fit_134976<-local({
    #Extract previously added sample data
    sampleInfo <-pData(gse_134976) 
    sampleInfo <- sampleInfo[,43:46]
    sampleInfo$dietcombo<-as.factor(rep(c("CC", "CO", "OC", "OO"), c(3,2,3,3)))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC", "CO", "OC", "OO")         #Rename Col names to make it prettier
    
    #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    contrasts <- makeContrasts(CC - OC, CO-OO, levels=design)
    
    ##Filter low expressed genes. Step omitted to keep all genes in analysis
    gse<-gse_134976
    # cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse_134976) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # #table(keep)                             # check how many genes are removed / retained.
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    #Use above specified function to produce fit
    fit <- lmFit(exprs(gse), design)
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit2 <- eBayes(fit2)
    fit_134976<-fit2
    
    anno <- fData(gse_134976)                    #Obtain a list with information for each probe
    anno <- dplyr::select(anno,GENENAME)    #Select columns that are of interest
    
    fit_134976$genes <- anno                    #Insert this new information into the table  with the results
    fit_134976
  })
  
  #Extract results with topTable and clean up omitting NA only selecting one sample per gene
  results_134976_CCvsOC<- local({
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results_F1ND <- topTable(fit_134976, coef=1, number=Inf, confint=T)%>%   #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")%>%                             #ProbeID to a column rather than Row name
      dplyr::rename(symbol=GENENAME)
    
    #Clean up full result sheet
    full_results_F1ND<-cleanup_results(full_results_F1ND)
  })
  
  results_134976_COvsOO<- local({
    #Creating a table with all genes by specifying the number argument for topTable
    full_results_F1HFD <- topTable(fit_134976, coef=2, number=Inf, confint=T)%>%   #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")%>%                             #ProbeID to a column rather than Row name
      dplyr::rename(symbol=GENENAME)
    
    #Clean up full result sheet
    full_results_F1HFD<-cleanup_results(full_results_F1HFD)
  })
  
  
  #Create a list with the TWO new result files
  GSE134976<-list(GSE134976_CCvsOC=results_134976_CCvsOC, 
                  GSE134976_COvsOO=results_134976_COvsOO)
}

#**********************************

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#GSE44901_CCvsOC&GSE44901_OOvsCO
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

GSE44901<-function(){
  
  #Download expression dataset
  gse_44901 <- getGEO(filename = "data/GSE44901_series_matrix.txt.gz")
  
  #Create fit including both conditions (F1HFD and F1ND)
  fit_44901<-local({
    
    gse<-gse_44901
    
    #exprs get the expression levels as a data frame and get the distribution
    #summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    ##Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             
    
    ## Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- sampleInfo[,1:2]              #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-sampleInfo%>%mutate(diet=stringr::str_extract(title, "^.{2}"))%>%
      transmute(dietcombo=case_when(diet=="LL"~"CC", diet=="WL"~"OC", diet=="WW"~"OO",
                                    diet=="LW"~"CO",TRUE~"else"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC", "OC", "CO", "OO")         #Rename Col names to make it prettier
    
    #specify contrasts
    contrasts <- makeContrasts(CC - OC, CO-OO, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    
    fit <- lmFit(exprs(gse), design)
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit_44901 <- eBayes(fit2)
    
    anno <- fData(gse)                    #Obtain a list with information for each probe
    anno <- dplyr::select(anno,Symbol)    #Select columns that are of interest
    
    fit_44901$genes <- anno                    #Insert this new information into the table  with the results
    fit_44901
  })
  
  #Extract results with topTable and clean up omitting NA only selecting one sample per gene
  results_44901_CCvsOC<- local({
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results_F1ND <- topTable(fit_44901, coef=1, number=Inf, confint=T)%>%   #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")                                          #ProbeID to a column rather than Row name
    
    #Clean up full result sheet
    full_results_clean_F1ND<-full_results_F1ND%>%
      dplyr::rename(symbol=Symbol)%>% #rename "Symbol" to match other datasets
      cleanup_results() 
    full_results_clean_F1ND
  })
  
  results_44901_COvsOO<- local({
    full_results_F1HFD <- topTable(fit_44901, coef=2, number=Inf, confint=T)%>%    #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")
    
    #Clean up full result sheet
    full_results_clean_F1HFD<-full_results_F1HFD%>%
      dplyr::rename(symbol=Symbol)%>%  #rename "Symbol" to match other datasets
      cleanup_results()  #select top signal probe for each gene
    full_results_clean_F1HFD
    
  })
  
  #Create a list with the TWO new result files
  GSE44901<-list(GSE44901_CCvsOC=results_44901_CCvsOC, 
                  GSE44901_COvsOO=results_44901_COvsOO)
}

#**********************************

#*#*#*#*#*#*#*#*#*#*#*
#GSE46359_M&GSE46359_F
#*#*#*#*#*#*#*#*#*#*#*

GSE46359<-function(){
  
  #Download expression dataset
  gse_46359 <- getGEO(filename = "data/GSE46359_series_matrix.txt.gz")
  
  #Create fit including both conditions (F1HFD and F1ND)
  fit_46359<-local({
    
    gse<-gse_46359
    
    #exprs get the expression levels as a data frame and get the distribution
    #summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    #exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    ##Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             
    
    ## Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- sampleInfo[,c(1,2,39)]              #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-sampleInfo%>%mutate(diet=stringr::str_extract(title,"[A-Z]{2,}"))%>%
      dplyr::rename(sex="gender:ch1")%>%
      transmute(dietcombo=case_when(diet=="LFD"&sex=="male"~"CC_M", diet=="WSD"&sex=="male"~"OC_M",
                                    diet=="LFD"&sex=="female"~"CC_F",diet=="WSD"&sex=="female"~"OC_F",
                                    TRUE~"else"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC_M", "OC_M", "CC_F", "OC_F")         #Rename Col names to make it prettier
    
    #specify contrasts
    contrasts <- makeContrasts(CC_M - OC_M, CC_F-OC_F, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    
    fit <- lmFit(exprs(gse), design)
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit_46359 <- eBayes(fit2)
    fit_46359
  })
  
  #Obtain annotation for samples
  anno_46359<-local({
    #Determine annotation of samples with the data stored under feature Data
    anno <- fData(gse_46359)%>%                     #Obtain a list with information for each probe
      dplyr::select(ID,gene_assignment)       #Select columns that are of interest
    
    #Select ENS items from the column gene information
    anno <- str_split(anno$gene_assignment, boundary("word"))%>%          #Split up the column in different strings/words
      enframe()%>%                                                        #Recreate a dataframe with a column-list
      cbind(anno)%>%                                                      #Combine with the original df to not use probe ID
      unnest(value)%>%dplyr::select(value, ID)%>%                         #Create a long format table storing all elements from "gene_id"
      dplyr::filter(str_detect(value, "^ENS")) 
    
    #match ENSMUSTXID to the corresponding gene symbols, command "columns(EnsDb.Mmusculus.v79)" shows available colnames
    gene_names<-select(EnsDb.Mmusculus.v79, keys(EnsDb.Mmusculus.v79), c("TXID", "GENEID", "GENENAME"))
    
    anno<-anno%>%left_join(gene_names, by=c("value"="TXID"))%>%dplyr::select(2,4)%>%
      na.omit(cols=c("GENENAME"))%>%  #filter out ID's where no symbol has been assigned
      mutate(ID = as.character(ID))%>%    #Pipe requires the mutate function, cannot use as.*** function in pipe
      distinct()%>% #only keep rows with a unique combination of ID and symbol
      #Select one symbol for each ID, no applied criteria for selection
      group_by(ID)%>%mutate(count=length(ID))%>% #group by ID and create a variable indicating how many symbols per ID
      slice_head(n=1)%>% #Select one symbol for each ID
      dplyr::select(1,2)%>% #Omit the previously created count variable
      dplyr::rename(symbol=GENENAME)
    anno
    
  })
  
  ?rename
  #Extract results with topTable and clean up omitting NA only selecting one sample per gene
  results_46359_M<- local({
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results_M <- topTable(fit_46359, coef=1, number=Inf, confint=T)%>%   #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")%>%                                         #ProbeID to a column rather than Row name
      left_join(anno_46359, by="ID")%>%   #Merging results with annotation previously created
      cleanup_results()   #Clean up full result sheet with previously created function
    
    full_results_M
  })
  
  results_46359_F<- local({
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results_F <- topTable(fit_46359, coef=2, number=Inf, confint=T)%>%   #Specify which contrast you want to look at (coef1=CC-OC, coef2=OC-OO)
      rownames_to_column("ID")%>%                                         #ProbeID to a column rather than Row name
      left_join(anno_46359, by="ID")%>%   #Merging results with annotation previously created
      cleanup_results()   #Clean up full result sheet with previously created function
    
    full_results_F
  })
  
  #Create a list with the TWO new result files
  GSE46359<-list(GSE46359_F=results_46359_F,
                 GSE46359_M=results_46359_M)
}

#**********************************

#*#*#*#*#
#GSE62715
#*#*#*#*#

GSE62715<-function(){
  
  #Import and transform expression data into an ExpressionSet
  gse_62715<-local({
    #Import file containing some sample and analysis information but no expression values
    gse_62715_pdata <- getGEO(filename = "data/GSE62715_series_matrix.txt.gz")
    
    #Import matrix containing expression data. Per above description the file contains "CPM and TMM" data
    #The counts table was downloaded from GREIN network (http://www.ilincs.org/apps/grein/?gse=GSE62715)
    GSE62715_cpm <- read.csv("data/GSE62715_GeneLevel_Normalized(CPM.and.TMM).csv.gz", row.names = 1)%>%
      dplyr::select(2:13) #select rows to be compatible with matrix, deselect symbols
    
    #Construct the Biobase ExpressionSet with the given elements
    gse_62715<-ExpressionSet(data.matrix(GSE62715_cpm)) #enter expression data in matrix form into the GSE
    pData(gse_62715)<-pData(gse_62715_pdata) #enter pData into the newly constructed expressionset
    
    #Enter fData with sample information as needed
    fData<- read.csv("data/GSE62715_GeneLevel_Normalized(CPM.and.TMM).csv.gz", row.names = 1)%>%
      dplyr::select(1)%>%   #only select column with the corresponding
      dplyr::rename(symbol=gene_symbol) #rename the var to correspond to other datasets later on
    
    fData(gse_62715)<-fData            #Enter the fData into the eSet
    
    #summary(exprs(gse_62715))    #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    exprs(gse_62715) <- log2(exprs(gse_62715))  #log2 transformation of the expression values
    #boxplot(exprs(gse_62715),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    gse_62715
  })
  
  #Produce fit
  fit_62715<-local({
    #Extract previously added sample data
    sampleInfo <-pData(gse_62715) 
    sampleInfo <- sampleInfo[,c(43,45)]
    sampleInfo$dietcombo<-as.factor(rep(c("CC", "OC"), c(6,6)))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$dietcombo)
    colnames(design) <- c("CC", "OC")         #Rename Col names to make it prettier
    
    #Specify the contrast of interest (the two differing maternal diets, same F1 diet)
    contrasts <- makeContrasts(CC - OC, levels=design)
    
    ##Selection of genes based on threshold criteria
    # gse<-gse_134976
    # cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse_134976) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # #table(keep)                             # check how many genes are removed / retained.
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    #Produce fit with the limma package, all genes included, no pre-selection
    fit <- lmFit(exprs(gse_62715), design)
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit_62715 <- eBayes(fit2)
    
    anno <- fData(gse_62715)%>%                    #Obtain a list with information for each probe
      dplyr::select(symbol)    #Select columns that are of interest
    
    fit_62715$genes <- anno                    #Insert this new information into the table  with the results
    fit_62715
  })
  
  #Extract results with topTable and clean up omitting NA only selecting one sample per gene
  results_62715<- local({
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results <- topTable(fit_62715, number=Inf, confint=T)%>%  
      rownames_to_column("ID")                    #ProbeID to a column rather than Row name
    
    #Clean up full result sheet
    results_62715<-cleanup_results(full_results)
  })
  
  #Create a list with the TWO new result files
  GSE62715<-list(GSE62715=results_62715)
}

#**********************************

#Create list with all the above functions####
list_results<-c(GSE133767(), GSE40903(), GSE123009_12w(), GSE123009_28w(),
                GSE134976(), GSE44901(), GSE46359(), GSE62715())

#Create Dataframe####
#Create a dataframe containing all above created datasets. Add a column containing the name of the GEO

df_results<-bind_rows(list_results, .id = "GEOSET")
write.csv(df_results, file="output/df_results.csv", row.names=FALSE)
