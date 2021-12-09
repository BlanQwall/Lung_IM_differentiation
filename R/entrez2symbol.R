entrez2symbol <- function(entrez_id, organism="hsa") {
  require(clusterProfiler)
  geneSymbol <- character(0)
  
  if (organism=="hsa") {
    orgdb <- "org.Hs.eg.db"
  } else {
    if (organism=="mmu") {
      orgdb <- "org.Mm.eg.db"
    }
  }
  
  for (i in 1:length(entrez_id)) {
    geneSymbol[i]  <- unlist(lapply(entrez_id[i], function(x) {
      tbl <- bitr(geneID=unlist(strsplit(x, split = "/")), 
                  fromType="ENTREZID", 
                  toType="SYMBOL", 
                  OrgDb = orgdb) 
      paste0(tbl$SYMBOL, collapse = "/")
    } ))
  }
  return(geneSymbol)
}



