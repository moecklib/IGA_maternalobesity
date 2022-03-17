#PLOTTING OF RAW EXPRESSION VALUES OF GSE164760
#Study the expression level of "immune" genes in NASH livers with and without
#cancer

#Beat Moeckli, University of Geneva, March 2022


#loading necessary packages
lapply(c("GEOquery", "ggprism", "ggbeeswarm","RColorBrewer", "dplyr", "tidyr", 
         "tidyverse", "limma", "EnsDb.Mmusculus.v79"), require, character.only = TRUE)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#Plot raw expression values for GSE164760####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#Define genes of interest
GenesOfInterest<-c("CXCL16", "CD4")

#Download expression dataset
gse_164760 <- getGEO(GEO="GSE164760")[[1]]

#Select columns that are of interest
anno<-fData(gse_164760)[, c("ID","Gene Symbol")]   #Select columns that are of interest
colnames(anno)<-c("ID","Gene_Symbol")

#Filter Affymetrix ID that match the genes of interest
affy_GOI<-anno[anno$Gene_Symbol%in%GenesOfInterest,1]

colnames(pData(gse_164760))
sample_info_pre<-data.frame(pData(gse_164760))
sample_info<-sample_info_pre[,"tissue.ch1",drop=F]

test<-pData(gse_164760)

#Filter the expression data according to the above defined genes of interest
Expr_GOI<-data.frame(exprs(gse_164760)[affy_GOI,])%>% #Extract the genes of interest
  tibble::rownames_to_column("ID")%>%         #Add row names to column
  left_join(anno, by="ID")%>%                 #Join with annotation file for gene symbol
  rowwise(c(ID, Gene_Symbol))%>%        #Group row_wise
  mutate(mean_expr=mean(c_across(where(is.numeric))))%>%  #Add mean expression level per proble
  relocate(c(ID, Gene_Symbol, mean_expr))%>%  #Rearrange columns
  ungroup()%>%                                
  arrange(desc(mean_expr))%>%                 #Arrange in descending order
  group_by(Gene_Symbol)%>%                    #Group by Gene symbol and select 
  slice_head(n=1)                   #the probe with the highest value

#Transpose columns
final_data_pre<-t(Expr_GOI)   #Transpose column
colnames(final_data_pre) <- final_data_pre[2,]  #Make second row col names
final_data_pre <- final_data_pre[-c(1,2), ]     #Remove first and second row

#Add sample characteristics and pivot to longer
final_data<-merge(final_data_pre, sample_info, by="row.names", all=T)%>%
  relocate(c(Row.names, tissue.ch1))%>%     #Change order of columns
  mutate(across(!c(Row.names, tissue.ch1), as.numeric))%>%
  pivot_longer(cols = where(is.numeric), names_to = "Gene", values_to = "Expr")

#Change name of the tissue type
final_data$tissue.ch1[final_data$tissue.ch1%in%"Non-tumoral NASH liver adjacent to HCC"]<-
  "NASH adjacent to HCC"
  
#Create plot with facet wrap
final_data%>%drop_na(tissue.ch1)%>%
ggplot(aes(x=tissue.ch1, y=Expr, fill=tissue.ch1))+
  geom_boxplot()+
  geom_beeswarm(shape=17)+
  facet_wrap(vars(Gene), nrow=1)+
  theme_prism()+
  theme(axis.text.x = element_text(angle = 45),
        axis.title.x = element_blank())

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Differential expression for GSE113508####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Analysis of a dataset comparing tumors in a steatotic environement and a non-
#steatotic hepatic environment.


#Download expression dataset
gse_113508 <- getGEO(GEO="GSE113508")[[1]]
  
fit_113508<-local({
    #Select the liver samples, the sample ID's have been identified with the GEO2R tool
    vehicle_samples<-c(c(1:3),c(7:9))        #Create a vector that corresponds to the liver and 12w samples
    gse <- gse_113508[ , vehicle_samples]       #selects only the liver samples
    
    #No log transformation needed as the distribution of the samples is between 2 and 14
    #exprs get the expression levels as a data frame and get the distribution
    #summary(exprs(gse))                 #Checks the distribution of the samples (summary), if larger than 16. a log2 transformation is needed
    
    #exprs(gse) <- log2(exprs(gse))      #log2 transformation of the expression values (one row per gene, one col per sample)
    #boxplot(exprs(gse),outline=FALSE)   #Plot the log transformed values, should be similar for all samples
    
    #Extraction of sample information which is stored under pData, this will help in determining which cols are helpful
    sampleInfo <-pData(gse)             

    # Select the columns which contain useful information about the samples for later comparison, rename for convenience
    sampleInfo <- dplyr::select(sampleInfo, "tumor growth condition:ch1")%>%
      dplyr::rename(liver_bckgrnd="tumor growth condition:ch1")              #Rename samples for convenience
    
    ##Prepare differential expression analysis with incorporating the levels to compare
    #Create a 4-level factor with incorporating the maternal and offspring diet in a single column
    #C=control diet, O=Obese/High Fat diet, the first letter stands for the F0 generation, the second for the F1 generation
    sampleInfo<-transmute(sampleInfo, liver_bckgrnd=case_when(liver_bckgrnd=="in healthy liver"~"healthy",
                                                              liver_bckgrnd=="in fatty liver"~"fatty",
                                                          TRUE~"NA"))
    
    #Differential expression, establish a model matrix with the 4-levels created previously (CC, OC, CO, OO)
    design <- model.matrix(~0+sampleInfo$liver_bckgrnd)
    colnames(design) <- c("fatty", "healthy")         #Rename Col names to make it prettier
    
    ##Filter low expressed genes from analysis
    # cutoff <- quantile(exprs(gse), 0.2)     # calculate 20th percentile (or median or else) expression level
    # is_expressed <- exprs(gse) > cutoff     # TRUE or FALSE for whether each gene is "expressed" in each sample
    # keep <- rowSums(is_expressed) > 2       # Identify genes expressed in more than 2 samples
    # gse <- gse[keep,]                       # subset to just those expressed genes, according to above criteria
    
    ##Determine the expression level in each pre-specified group
    fit <- lmFit(exprs(gse), design)
    
    #specify contrasts
    contrasts <- makeContrasts(healthy - fatty, levels=design)     #Specify the two contrasts of interest (the two differing maternal diets, same F1 diet)
    
    fit2 <- contrasts.fit(fit, contrasts)           #Apply empirical Bayes step to obtain DE statistics and P-values
    fit2 <- eBayes(fit2)
    fit_113508<-fit2
  })
  

  #Create output with annotation
  results_113508<-local({
    
    anno<-local({
    ##Extract the annotation data from the affymetrix description sheet
    #Determine annotation of samples with the data stored under feature Data
    anno <- fData(gse_113508)%>%                     #Obtain a list with information for each probe
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
    anno<-anno%>%na.omit(symbol)
    anno
    })
    
    #Insert this new information into the table  with the results
    fit_113508$genes <- anno
    
    #Creating a table with all genes by specifying the number argument for topTable
    full_results<- topTable(fit_113508, number=Inf, confint = T)%>%    #Specify which contrast you want to look if several
      na.omit(ID)              #ProbeID to a column rather than Row name
    
    #Clean up full result sheet
    full_results<-full_results%>%
      na.omit(cols=c("symbol", "logFC"))%>%  #omit rows/probes with no symbol assigned
      group_by(symbol)%>%
      arrange(desc(AveExpr))%>%
      slice_head(n = 1)  #select top signal probe for each gene
    full_results
  })
  
  write.csv(results_113508, file="LimmaResults_GSE11358.csv")
  