entrez2symbol <- function(entrez_id) {
  require(clusterProfiler)
  geneSymbol <- character(0)
  for (i in 1:length(entrez_id)) {
    geneSymbol[i]  <- unlist(lapply(entrez_id[i], function(x) {
      tbl <- bitr(geneID=unlist(strsplit(x, split = "/")), fromType="ENTREZID", toType="SYMBOL", OrgDb = "org.Hs.eg.db") 
      paste0(tbl$SYMBOL, collapse = "/")
    } ))
  }
  return(geneSymbol)
}



