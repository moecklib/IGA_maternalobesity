


source("code/pathway_analysis/lib_reactome.R")


# Load reactome annotation
R <- reactome("Mus musculus")

# Output Reactome annotations for Beat
local({
  A <- merge(R$mapping_lowest[c("reactome_id","entrez_id")],R$ancestors,by.x="reactome_id",by.y="descendant_id")
  A <- A[c("entrez_id","ancestor_id")]
  A <- unique(A)
  A$symbol <- select(org.Mm.eg.db,sub("^G:","",A$entrez_id),columns = "SYMBOL")$SYMBOL
  stopifnot(all(!duplicated(V(R$hierachy)$pathway))) # Check pathway names are unique
  A$pathway <- V(R$hierachy)[A$ancestor_id]$pathway
  PW <- splitAsList(paste0(sub("^R:","",A$ancestor_id),":",A$pathway),A$symbol)
  PW <- stack(PW)
  PW <- subset(PW,!value %in% "root:root")
  write.csv(PW,row.names=FALSE,file=gzfile("output/reactome_pathways.csv.gz"))
})