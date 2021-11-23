replaceEntrezID <- function(enrichResults) {
  result <- enrichResults@result
  result$geneID <- entrez2symbol(result$geneID)
  enrichResults@result <- result
  return(enrichResults)
}