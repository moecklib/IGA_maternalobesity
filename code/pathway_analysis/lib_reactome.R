


library(igraph)
library(org.Mm.eg.db)



#' load an igraph containing metadata of the pathways and the hierarchy of the reactome pathways
reactome_hierachy <- function(org="Homo sapiens") {
  x <- read.table("data/reactome/ReactomePathways.txt.gz",sep="\t",quote="",col.names = c("id","pathway","org"),comment="")
  h <- read.table("data/reactome/ReactomePathwaysRelation.txt.gz",sep="\t",quote="",col.names = c("parent","child"),comment="")
  x$id <- paste0("R:",x$id)
  h$parent <- paste0("R:",h$parent)
  h$child <- paste0("R:",h$child)
  
  # check data consistency
  stopifnot(all(h$parent %in% x$id))
  stopifnot(all(h$child %in% x$id))
  stopifnot(all(x$org[match(h$child,x$id)]==x$org[match(h$parent,x$id)]))
  
  # subset on selected organism
  x <- x[x$org %in% org,]
  h <- h[(h$parent %in% x$id) & (h$child %in% x$id),]
  
  # create the graph
  x$reactome_id <- x$id
  x$pathway <- ifelse(x$pathway %in% x$pathway[duplicated(x$pathway)],paste0(x$reactome_id,"|",x$pathway),x$pathway)
  g <- graph_from_data_frame(h,TRUE,x)
  
  # add a root node
  g <- g %>% 
    add_vertices(1L,attr=list(name="root",pathway="root")) %>%
    add_edges(rbind("root",V(g)[!.to(E(g))]$name))
  
  g
}

reactome_ancestors <- function(h) {
  lnk <- lapply(V(h),function(v) subcomponent(h,v,mode="in")$name)
  lnk <- stack(lnk)
  names(lnk) <- c(ind="descendant_id",values="ancestor_id")[names(lnk)]
  lnk
}

# load annotations
reactome_mapping_lowest <- function(h) {
  a <- read.table("data/reactome/NCBI2Reactome.txt.gz",sep="\t",quote="",comment="",col.names = c("entrez_id","reactome_id","url","event","evidence","org"))
  a <- a[c("reactome_id","entrez_id","event","evidence","org")]
  a$reactome_id <- paste0("R:",a$reactome_id)
  a$entrez_id <- paste0("G:",a$entrez_id)
  a <- a[a$reactome_id %in% V(h)$reactome_id,]
  stopifnot(all(V(h)[a$reactome_id]$org == a$org)) # check organism are compatible
  a
}


reactome <- function(...) {
  h <- reactome_hierachy(...)
  list(
    hierachy = h,
    mapping_lowest = reactome_mapping_lowest(h),
    ancestors = reactome_ancestors(h)
  )
}


reactome_all_levels <- function(R=reactome())  {
  A <- merge(R$mapping_lowest[c("reactome_id","entrez_id")],R$ancestors,by.x="reactome_id",by.y="descendant_id")
  A <- A[c("entrez_id","ancestor_id")]
  A <- unique(A)
  A$symbol <- select(org.Mm.eg.db,sub("^G:","",A$entrez_id),columns = "SYMBOL")$SYMBOL
  stopifnot(all(!duplicated(V(R$hierachy)$pathway))) # Check pathway names are unique
  A$pathway <- V(R$hierachy)[A$ancestor_id]$pathway
  splitAsList(A$pathway,A$symbol)
}


tree_from_dag <- function(dag,root) {
  stopifnot(is_dag(dag))
  leafs <- V(dag)[!.from(E(dag))]
  paths <- all_simple_paths(dag,root,to=leafs,mode = "out")
  
  n <- paths %>%
    lapply(Reduce,f=function(a,b) paste0(a,"/",b),accumulate=TRUE) %>%
    unlist() %>%
    unique()
  e <- cbind(dirname(n),n)
  e <- e[e[,1] != ".",]
  g <- graph_from_edgelist(e)
  revmap <- as.integer(basename(V(g)$name))
  V(g)$name <- paste0("P:",V(g)$name)
  
  list(tree = g,vertex_index = revmap)
}




reactome_hsa <- function(fc,root="root",depth=2,R=reactome())  {
  # Select reactome nodes according to parameters
  H <- make_ego_graph(R$hierachy,depth,root,mode = "out")[[1]]
  U <- tree_from_dag(H,root)
  V(U$tree)$pathway <- V(H)$pathway[U$vertex_index]
  V(U$tree)$reactome_id <- V(H)$reactome_id[U$vertex_index]
  #ggraph(U$tree, 'tree',circular=TRUE) + geom_edge_link() + geom_node_point() + coord_equal()
  
  # Remap annotations
  leafs <- unique(V(U$tree)[!.from(E(U$tree))]$reactome_id)
  leafs_remap <- R$ancestors[R$ancestors$ancestor_id %in% leafs,]
  leafs_remap <- merge(leafs_remap,R$mapping_lowest[c("reactome_id","entrez_id")],by.x="descendant_id",by.y="reactome_id")
  leafs_remap$reactome_id <- leafs_remap$ancestor_id
  A <- unique(rbind(
    R$mapping_lowest[R$mapping_lowest$reactome_id %in% V(U$tree)$reactome_id,c("reactome_id","entrez_id")],
    leafs_remap[c("reactome_id","entrez_id")]
  ))
  A <- merge(A,data.frame(tree_vid=as_ids(V(U$tree)),reactome_id=V(U$tree)$reactome_id))
  A <- A[!duplicated(A[c("tree_vid","entrez_id")]),]
  
  # Extract gene vertices
  G <- unique(data.frame(name=A$entrez_id,entrez_id=A$entrez_id))
  G$symbol <- select(org.Hs.eg.db,sub("^G:","",G$entrez_id),columns = "SYMBOL")$SYMBOL
  
  # Merge information into a single graph
  H2 <- U$tree %>%
    add_vertices(nrow(G),attr=as.data.frame(G)) %>%
    add_edges(rbind(A$tree_vid,A$entrez_id))
  
  # add fc and remove genes with missing values
  V(H2)$fc <- fc[match(V(H2)$symbol,names(fc))]
  H2 <- H2 - V(H2)[grepl("^G:",name) & is.na(fc)]
  
  H2
}




