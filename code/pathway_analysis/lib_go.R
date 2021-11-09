

library(igraph)
library(Matrix)
library(IRanges)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Gene Ontology
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# read OBO format and return an igraph with all "is_a" and "part_of" relations
read_obo <- function(obo_file) {
  pat <- "^([^:]+) *: *(.*)$"
  obo <- readLines(obo_file)
  obo <- obo[grepl("^\\[",obo) | grepl(pat,obo)]
  obo <- DataFrame(
    section = cumsum(grepl("^\\[",obo)),
    key = sub(pat,"\\1",obo),
    value = sub(pat,"\\2",obo)
  )
  obo <- subset(obo,section %in% section[key=="[Term]"])
  obo <- subset(obo,key %in% c("is_a","name","id","is_obsolete","namespace","relationship"))
  
  get_relationship <- function(obo) {
    R <- with(obo[obo$key=="relationship",],CharacterList(split(value,section)))
    R <- sub(" *!.*$","",R) # remove comments
    pat <- "^([^ ]+) (GO:[0-9]+)$"
    stopifnot(all(grepl(pat,unlist(R))))
    DataFrame(
      relationship_type = sub(pat,"\\1",R),
      relationship_goid = sub(pat,"\\2",R)
    )
  }
  
  obo$section <- factor(obo$section)
  obo <- DataFrame(
    goid = with(subset(obo,key=="id"),drop(CharacterList(split(value,section)))),
    term = with(subset(obo,key=="name"),drop(CharacterList(split(value,section)))),
    is_obsolete = with(subset(obo,key=="is_obsolete"),drop(CharacterList(split(value,section)))) %in% "true",
    namespace = with(subset(obo,key=="namespace"),drop(CharacterList(split(value,section)))),
    is_a = with(subset(obo,key=="is_a"),CharacterList(split(value,section))),
    get_relationship(obo)
  )
  obo$is_a <- sub(" *!.*$","",obo$is_a)
  
  is_a <- stack(setNames(obo$is_a,obo$goid))
  is_a$type <- "is_a"
  relationship <- stack(setNames(obo$relationship_goid,obo$goid))
  relationship$type <- unlist(obo$relationship_type)
  
  go <- graph_from_data_frame(rbind(is_a,relationship),vertices = as.data.frame(obo[1:4]))
  V(go)$term_and_id <- sprintf("%s: %s",V(go)$name,V(go)$term)
  go
}


query_go_basic <- function() {
  read_obo("http://purl.obolibrary.org/obo/go/go-basic.obo")
}


# Given a graph, compute the sparse ancestor binary matrix that associate to each node all its ancestor nodes (including itself)
# WARNING: I suspect the function to loop infinitely if the given graph contains cycles !
ancestors <- function(go) {
  stopifnot(is.dag(go))
  m0 <- as_adjacency_matrix(go,sparse=TRUE)>0
  diag(m0) <- TRUE
  m <- m0
  n <- 0
  while (n!=sum(m)) {
    n <- sum(m)
    m <- (m + (m %*% m0))>0
  }
  m  
}


goall <- function(go) {
  goall <- ancestors(go)
  ij <- Matrix::which(goall,arr.ind=TRUE)
  splitAsList(colnames(goall)[ij[,2]],rownames(goall)[ij[,1]])
}

# simplify a graph
igraph_remove_node_keeping_connections <- function(G,v) {
  e <- as.vector(t(as.matrix(expand.grid(V(G)[.innei(v)]$name,V(G)[.outnei(v)]$name))))
  add_edges(G-v,e)
}


# simplify the graph by removing transitive edges
igraph_remove_transitive_edges <- function(G) {
  E(G)$is_transitive <- NA
  for(i in seq_along(E(G))) {
    e <- E(G)[i]
    g <- G - e
    E(G)$is_transitive[i] <- is.finite(distances(g,tail_of(G,e)$name,head_of(G,e)$name,mode = "out"))
  }
  G - E(G)[is_transitive]
}




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Uniprot
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


# Retreive for Uniprot all GO annotations of the given organism
query_uniprot <- function(org) {
  url <- "https://www.uniprot.org/uniprot/?format=tab&columns=id,entry%20name,protein%20names,genes(PREFERRED),genes,ec,go-id,keyword-id"
  url <- paste0(url,"&query=organism:",org)
  uniprot <- read.table(url,sep="\t",comment="",quote="",stringsAsFactors=FALSE,header=TRUE,check.names=FALSE)
  uniprot <- as(uniprot,"DataFrame")
  
  uniprot$Gene_names <- CharacterList(strsplit(uniprot$"Gene names","( +)|( */ *)"))
  uniprot$Gene_ontology_IDs <- CharacterList(strsplit(uniprot$"Gene ontology IDs"," *; *"))
  uniprot$EC_number <- CharacterList(strsplit(uniprot$"EC number"," *; *"))
  uniprot$Keyword_IDs <- CharacterList(strsplit(uniprot$"Keyword ID"," *; *"))
  within(uniprot,`Gene names` <- `Gene ontology IDs` <- `EC number` <- `Keyword ID` <- NULL)
}



query_uniprot_keywords <- function() {
  lnk <- url('https://www.uniprot.org/keywords/?query=*&format=tab&force=true&columns=id&compress=yes')
  uniprot_keywords <- read.table(gzcon(lnk,text = TRUE),sep="\t",header=TRUE,comment="",check.names = FALSE,stringsAsFactors = FALSE,quote="")
  uniprot_keywords <- as(uniprot_keywords,"DataFrame")
  uniprot_keywords
}






#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# GSEA
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


#' Test GOterm enrichement for the given list of genes
#' @param query a character vector of gene names matching with names(genes_annotations)
#' @param genes_annotations a CharacterList that list for each gene the sets it belong to. 
#'        The list should be limited to the gene of the universe, for example genes expressed in the brain.
gene_set_enrichment_analysis <- function(query,genes_annotations) {
  #query <- x$locus_tag[x$cluster_id %in% 14]
  #genes_annotations <- x$goall_id
  
  # remove duplicated and only keep genes that are not annotated
  genes_annotations <- unique(genes_annotations)
  genes_annotations <- genes_annotations[lengths(genes_annotations)>0,]
  
  query <- intersect(query,names(genes_annotations))
  
  gsea <- with(stack(genes_annotations),splitAsList(as.character(name),value))
  gsea <- DataFrame(
    geneset_name = names(gsea),
    geneset_genes = gsea,
    geneset_query_genes = gsea[gsea %in% query]
  )
  gsea$geneset_size <- lengths(gsea$geneset_genes)
  gsea$n <- lengths(gsea$geneset_query_genes)
  gsea$query_size <- length(query)
  gsea$universe_size <- length(genes_annotations)
  gsea$expected_n <- gsea$query_size * gsea$geneset_size / gsea$universe_size
  gsea$fold_enrichment <- gsea$n / gsea$expected_n
  gsea$enrichment_pval <- 1 - phyper(gsea$n-1L,gsea$geneset_size,gsea$universe_size-gsea$geneset_size,gsea$query_size)
  gsea$enrichment_qval <- p.adjust(gsea$enrichment_pval,"fdr")
  gsea <- gsea[order(gsea$enrichment_pval),]
  gsea
}





#' Test GOterm enrichement for the given list of genes
#' @param query a logical matrix of selected genes (rownames are gene symbols, colnames are names of the geneset)
#' @param genes_annotations a CharacterList that list for each gene the sets it belong to. 
#'        The list should be limited to the gene of the universe, for example genes expressed in the brain.
gene_set_enrichment_analysis_v2 <- function(query,genes_annotations) {
  #query <- BS;genes_annotations <- reactome_all_levels()
  stopifnot(is.matrix(query))
  stopifnot(is.logical(query))
  
  # remove duplicates and only keep genes that are annotated
  genes_annotations <- unique(genes_annotations)
  genes_annotations <- genes_annotations[lengths(genes_annotations)>0,]
  genes_annotations <- genes_annotations[names(genes_annotations) %in% rownames(query)]
  
  query <- query[rownames(query) %in% names(genes_annotations),,drop=FALSE]
  
  gsea <- with(stack(genes_annotations),splitAsList(as.character(name),value))
  gsea <- DataFrame(
    geneset_name = names(gsea),
    geneset_genes = gsea
  )
  gsea$geneset_size <- lengths(gsea$geneset_genes)
  
  gsea$overlap_genes <- as(apply(query,2,function(v) {
    gsea$geneset_genes[gsea$geneset_genes %in% rownames(query)[v]]
  }),"DataFrame")
  gsea$overlap_size <- sapply(gsea$overlap_genes,lengths)
  gsea$query_size <- matrix(rep(colSums(query),each=nrow(gsea)),nrow(gsea))
  gsea$universe_size <- matrix(rep(colSums(!is.na(query)),each=nrow(gsea)),nrow(gsea))
  
  gsea$expected_overlap_size <- gsea$query_size * gsea$geneset_size / gsea$universe_size
  gsea$fold_enrichment <- gsea$overlap_size / gsea$expected_overlap_size
  gsea$enrichment_pval <- 1 - phyper(gsea$overlap_size-1L,gsea$geneset_size,gsea$universe_size-gsea$geneset_size,gsea$query_size)
  gsea$enrichment_qval <- array(p.adjust(gsea$enrichment_pval,"fdr"),dim(gsea$enrichment_pval),dimnames(gsea$enrichment_pval)) 
  gsea <- gsea[order(rowMin(gsea$enrichment_pval)),]
  gsea
}







#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Plotting and selection methods
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

gsea_df <- function(gene_sets,locus_tags,cluster_ids) {
  cluster_ids <- factor(cluster_ids)
  gene_sets <- gene_sets[names(gene_sets) %in% locus_tags[!is.na(cluster_ids)]]
  gsea <- SplitDataFrameList(lapply(levels(cluster_ids),function(cluster_id) {
    gsea <- gene_set_enrichment_analysis(
      locus_tags[cluster_ids %in% cluster_id],
      gene_sets
    )
    gsea$rank <- seq_along(gsea$geneset_name)
    gsea
  }))
  names(gsea) <- levels(cluster_ids)
  gsea <- stack(gsea,"cluster_id")
  gsea
}


select_topset <- function(GSEA,pval=1e-4) {
  GSEA$geneset_name <- paste0(GSEA$geneset_name,"/",GSEA$geneset_size)
  if (pval < 1) {
    GSEA <- GSEA[GSEA$geneset_name %in% GSEA$geneset_name[GSEA$enrichment_pval <= pval],]
    GSEA$geneset_name <- factor(GSEA$geneset_name,with(GSEA[GSEA$enrichment_pval <= pval,],unique(geneset_name[order(cluster_id,rank)])))
  } else {
    GSEA <- GSEA[GSEA$geneset_name %in% GSEA$geneset_name[GSEA$rank <= pval],]
    GSEA$geneset_name <- factor(GSEA$geneset_name,with(GSEA[GSEA$rank <= pval,],unique(geneset_name[order(cluster_id,rank)])))
  }
  GSEA$phred <- pmin(-log10(GSEA$enrichment_pval),5)
  GSEA
}

gg_gsea_area <- function(GSEA) {
  ggplot(as.data.frame(GSEA[c("geneset_name","enrichment_pval","rank","n","geneset_size","cluster_id","phred")])) +
    geom_point(aes(x=cluster_id,y=geneset_name,size=phred,color=cluster_id)) +
    scale_size_area() +
    theme_minimal() +
    theme(panel.grid = element_blank()) + xlab("cluster") + ylab("") +
    guides(color=guide_none())
}

gg_gsea_bitmap <- function(GSEA) {
  X <- stack(setNames(GSEA$geneset_genes,GSEA$geneset_name))
  X$cluster_order <- mcols(rna)$cluster_order[match(X$value,rownames(rna))]
  X$cluster_id <- mcols(rna)$cluster_id[match(X$value,rownames(rna))]
  ggplot(as(X,"data.frame")) +
    geom_point(aes(x=name,y=cluster_order,color=factor(cluster_id)),size=0.1) +
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
}




