replaceEntrezID <- function(enrichResults, organism="hsa") {
  result <- enrichResults@result
  result$geneID <- entrez2symbol(result$geneID, organism = organism)
  enrichResults@result <- result
  return(enrichResults)
}