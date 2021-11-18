
library(ggplot2)
library(patchwork)
library(ggraph)
library(SummarizedExperiment)
source("code/pathway_analysis/lib_gsea.R")
source("code/pathway_analysis/lib_obo.R")


# heatmap
gg_avg_heatmap <- function(avg) {
  avg[is.na(avg)] <- 0
  avg[avg<=-1] <- -1
  avg[avg>=+1] <- +1
  avg <- avg[hclust(dist(avg),method="ward.D")$order,]
  z <- data.frame(
    symbol = factor(rownames(avg),rownames(avg))[row(avg)],
    GEOSET = factor(colnames(avg),colnames(avg))[col(avg)],
    logFC=as.vector(avg)
  )
  ggplot(z) +
    geom_tile(aes(x=GEOSET,y=symbol,fill=logFC)) +
    scale_fill_gradient2(low="#C900C9",high="yellow",mid="black",midpoint=0) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5),
          legend.position = "top"
    )
}

# Load data compiled by Beat
x <- local({
  x <- read.csv("output/df_results.csv",check.names = FALSE)
  stopifnot(max(table(x$symbol,x$GEOSET)) == 1)
  SummarizedExperiment(
    lapply(x[-1:-3],function(v) tapply(v,list(x$symbol,x$GEOSET),mean,na.rm=TRUE))
  )
})

# Add reactome annotations
mcols(x)$reactome <- local({
  A <- read.csv("output/reactome_pathways.csv.gz")
  splitAsList(A$value,factor(A$name,rownames(x)))
})

mcols(x)$msigdb <- local({
  A <- read_gmt("data/msigdb.v7.4.symbols.gmt.gz")
  with(stack(setNames(A$MEMBERS,A$STANDARD_NAME)),splitAsList(as.character(name),factor(value,make.unique(toupper(rownames(x))))))
})


# Filter data
x <- x[,setdiff(colnames(x),c("GSE133767","GSE62715"))]
x <- x[rowSums(!is.na(assay(x))) > 0]

# Show an heatmap of genes seen in almost all dataset but 2
gg_avg_heatmap(assay(x[rowSums(is.na(assay(x))) <= 2]))
ggsave("output/heatmap.png",width = 2,height=10)
ggsave("output/heatmap.pdf",width = 2,height=10)




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Test Msigdb
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
X <- x[lengths(mcols(x)$msigdb)>0]

# Select annotated genes with low P-value
set.seed(123)
mcols(X)$ud <- data.frame(
  dn = local({
    P <- -log10(assay(X,"P.Value"))
    P <- P * ifelse(assay(X,"logFC")<0,-1,1)
    R <- apply(P,2,rank,ties.method="random")
    R[is.na(P)] <- NA
    R <- (R <= 250)
    R & (assay(X,"logFC") <= -0.1)
  }),
  up = local({
    P <- -log10(assay(X,"P.Value"))
    P <- P * ifelse(assay(X,"logFC")<0,-1,1)
    R <- apply(-P,2,rank,ties.method="random")
    R[is.na(P)] <- NA
    R <- (R <= 250)
    R & (assay(X,"logFC") >= +0.1)  })
)
colSums(mcols(X)$ud,na.rm=TRUE)
colSums(!is.na(mcols(X)$ud),na.rm=TRUE)
gsea <- multi_gsea(as.matrix(mcols(X)$ud),mcols(X)$msigdb)

GSEA <- gsea[1:30,]



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Now Reactome
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Select annotated genes with low P-value
X <- x[lengths(mcols(x)$reactome)>0]
set.seed(123)
mcols(X)$ud <- data.frame(
  dn = local({
    P <- -log10(assay(X,"P.Value"))
    P <- P * ifelse(assay(X,"logFC")<0,-1,1)
    R <- apply(P,2,rank,ties.method="random")
    R[is.na(P)] <- NA
    R <- (R <= 250)
    R & (assay(X,"logFC") <= -0.1)
  }),
  up = local({
    P <- -log10(assay(X,"P.Value"))
    P <- P * ifelse(assay(X,"logFC")<0,-1,1)
    R <- apply(-P,2,rank,ties.method="random")
    R[is.na(P)] <- NA
    R <- (R <= 250)
    R & (assay(X,"logFC") >= +0.1)  })
)
gg_avg_heatmap(assay(X[rowSums(mcols(X)$ud,na.rm=TRUE)>0,]))
ggsave("output/heatmap_top250.png",width = 1.5,height=5)
ggsave("output/heatmap_top250.pdf",width = 1.5,height=5)


# Perform enrichment analysis on all sets
gsea <- multi_gsea(as.matrix(mcols(X)$ud),mcols(X)$reactome)
gsea$reactome_id <- sub(":.*","",gsea$geneset_name)
gsea$reactome_id <- paste0("R:",gsea$reactome_id)





# Select enriched sets
GSEA <- gsea[rowSums(gsea$enrichment_pval <= 1e-2)>=3,]


source("code/pathway_analysis/lib_reactome.R")
R <- reactome("Mus musculus")
local({
  GSEA$reactome_ancestors_id <- splitAsList(R$ancestors$ancestor_id,R$ancestors$descendant_id)[GSEA$reactome_id]
  GSEA$lev1 <- local({
    lev1 <- ego(R$hierachy,1,"root")[[1]]
    GSEA$reactome_ancestors_id[GSEA$reactome_ancestors_id %in% setdiff(lev1$name,"root")]  
  })
  #GSEA <- GSEA[!GSEA$reactome_id %in% unlist(GSEA$lev1),]
  GSEA$lev1 <- relist(V(R$hierachy)[unlist(GSEA$lev1)]$pathway,GSEA$lev1)
  GSEA$lev1 <- unique(GSEA$lev1)
  GSEA$lev1 <- unstrsplit(GSEA$lev1," | ")
  
  M <- GSEA$enrichment_pval
  M <- reshape2::melt(M)
  M$pathway <- GSEA$pathway[M$Var1]
  M$lev1 <- GSEA$lev1[M$Var1]
  M$sign <- ifelse(sub("^(up|dn)\\.(.*)$","\\1",M$Var2)%in%"dn",-1,+1)
  M$grp <- sub("^(up|dn)\\.(.*)$","\\2",M$Var2)
  ggplot(M) + 
    facet_grid(lev1~grp,scales = "free_y",space = "free_y") +
    geom_col(aes(x=Var1,y=sign * -log10(pmax(value,1e-6)),fill=factor(sign)),color="black") +
    coord_flip() +
    scale_fill_manual(values = c("#990099","yellow")) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),legend.position = "none") +
    ylab("phred") + xlab("")
})
ggsave("output/barplot_all.pdf",width=18,height=8)





source("code/pathway_analysis/lib_util.R")
library(ggraph)

d <- local({
  # Full reactome dendrogram
  d <- tree_from_dag(R$hierachy,"root")
  V(d$tree)$reactome_id <- V(R$hierachy)$reactome_id[d$vertex_index]
  d <- d$tree
  V(d)$pathway <- V(R$hierachy)$pathway[match(V(d)$reactome_id,V(R$hierachy)$reactome_id)]
  V(d)$depth <- nchar(gsub("[^/]","",V(d)$name))

  # reduce to GSEA sets only
  d <- induced_subgraph(d,na.omit(bfs(d,V(d)[reactome_id %in% GSEA$reactome_id],neimode = "in",unreachable=FALSE)$order))
  
  # Simplify the tree (remove intermediate node with single links)
  while(TRUE) {
    singletons <- V(d)[(degree(d,mode = "in") == 1L) & (degree(d,mode = "out") == 1L)]
    if (length(singletons)==0) break
    d <- d %>%
      add_edges(c(V(d)[.innei(singletons[1])],V(d)[.outnei(singletons[1])])) %>%
      delete_vertices(singletons[1]$name)
  }
  d
})
V(d)$terminal <- FALSE;V(d)[!.from(E(d))]$terminal <- TRUE





# Show simplified dendrogram
p0 <- ggraph(d,layout="dendrogram") +
  geom_edge_bend() +
  geom_node_label(aes(label=str_break(pathway,20),filter=!leaf),size=3,hjust=0.5,vjust=0.5) +
  #geom_node_text(aes(label=pathway,filter=leaf),size=3,hjust=0,vjust=0.5) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  scale_x_continuous(expand=c(0,1.5,0,1.5)) +
  coord_flip()



p1 <- local({
  # Get enrichement p_values for all experiments
  vd_pval <- GSEA$enrichment_pval[match(V(d)$reactome_id,GSEA$reactome_id),]
  V(d)$x_layout <- p0$data$x[match(V(d)$name,p0$data$name)]
  M <- reshape2::melt(vd_pval)
  M$sign <- ifelse(sub("^(up|dn)\\.(.*)$","\\1",M$Var2)%in%"dn",-1,+1)
  M$grp <- sub("^(up|dn)\\.(.*)$","\\2",M$Var2)
  M <- merge(M,within(as_data_frame(d,"v"),gsea_id<-paste0(sub("^R:","",reactome_id),":",pathway)),by.x="Var1",by.y="gsea_id")  # Add layout from p0
  M <- M[M$terminal,]
  
  # coordinaltes of horizontal lines
  hl <- sapply(V(d)[.outnei(V(d)[1])]$name,function(v){vx <- subcomponent(d,v,mode = "out")$x_layout;range(vx+0.5,vx-0.5)})
  hl <- unique(as.vector(hl))
  
  ggplot(M) + facet_grid(.~grp) +
    geom_rect(aes(ymin=x_layout-0.4,ymax=x_layout+0.4,xmin=0,xmax=sign * -log10(pmax(value,1e-6)),fill=factor(sign)),color="black") +
    scale_fill_manual(values = c("#990099","yellow")) +
    scale_y_continuous(breaks=V(d)[terminal]$x_layout,labels=V(d)[terminal]$pathway,expand=c(0,1,0,1)) +
    geom_hline(yintercept=hl) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),legend.position = "none") +
    xlab("phred")
})

p1 + p0 + plot_layout(width=c(3,1))
ggsave("output/barplot2.pdf",width=17,height=5)



local({
  gsea$geneset_genes <- unstrsplit(gsea$geneset_genes,",")
  gsea$overlap_genes <- endoapply(gsea$overlap_genes,unstrsplit,sep=",")
  gsea$query_size <- gsea$universe_size <- gsea$enrichment_qval <- gsea$reactome_ancestors_id <- gsea$ fold_enrichment <- gsea$expected_overlap_size <- NULL
  write.table(as.data.frame(gsea),sep="\t",file = "output/barplot_all.tsv",row.names = FALSE)
})






