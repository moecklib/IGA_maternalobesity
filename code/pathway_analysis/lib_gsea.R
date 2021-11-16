

library(IRanges)
library(ggplot2)

#' Test GOterm enrichement for the given list of genes
#' @param query a logical matrix of selected genes (rownames are gene symbols, colnames are names of the geneset)
#' @param genes_annotations a CharacterList that list for each gene the sets it belong to. 
#'        The list should be limited to the gene of the universe, for example genes expressed in the brain.
multi_gsea <- function(query,genes_annotations) {
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
    gsea$query_size <- matrix(rep(colSums(query,na.rm=TRUE),each=nrow(gsea)),nrow(gsea))
    gsea$universe_size <- matrix(rep(colSums(!is.na(query)),each=nrow(gsea)),nrow(gsea))
    
    gsea$expected_overlap_size <- gsea$query_size * gsea$geneset_size / gsea$universe_size
    gsea$fold_enrichment <- gsea$overlap_size / gsea$expected_overlap_size
    gsea$enrichment_pval <- 1 - phyper(gsea$overlap_size-1L,gsea$geneset_size,gsea$universe_size-gsea$geneset_size,gsea$query_size)
    gsea$enrichment_qval <- array(p.adjust(gsea$enrichment_pval,"fdr"),dim(gsea$enrichment_pval),dimnames(gsea$enrichment_pval)) 
    
    # sort the output
    o <- gsea$enrichment_pval
    o[is.na(o)] <- +Inf
    o <- order(Biobase::rowMin(o))
    gsea[o,]
}



read_gmt <- function(con) {
  txt <- readLines(con)
  pat <- "^([^\t]*)\t([^\t]*)\t(.*)$"
  z <- DataFrame(
    STANDARD_NAME = sub(pat,"\\1",txt),
    EXTERNAL_DETAILS_URL = sub(pat,"\\2",txt),
    MEMBERS = CharacterList(strsplit(sub(pat,"\\3",txt),"\t",fixed=TRUE))
  )
}


gg_matrix_hcol <- function(m,fill,grp) {
  stopifnot(ncol(m)==length(fill))
  stopifnot(length(grp)==length(fill))
            
  p <- reshape2::melt(m)
  p$fill <- fill[p$Var2]
  p$grp <- grp[p$Var2]
  ggplot(p) + facet_wrap(.~grp,nrow=1) +
    geom_col(aes(x=Var1,y=value,fill=fill)) + 
    coord_flip() +
    geom_hline(yintercept=0) +
    theme(panel.spacing = unit(0,"mm"),panel.grid = element_blank(),strip.text = element_text(angle=90,hjust=0)) + 
    xlab("") + ylab("")
}





