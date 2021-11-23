Seurat2CellFreqTable <- function(seuratObject, slotName=NULL, cells = NULL, decreasing = TRUE, top = FALSE, return.cells = FALSE) {
  
  require(Seurat)
  
  # load data: 
  data <- seuratObject@meta.data
  # if cells are defined, limit data: 
  if (!is.null(cells)) {data <- data[cells,]}
  if (is.null(slotName)) { data <- seuratObject@active.ident
    slotName <- "active.ident"
    
  } else {
  cellNames <- rownames(data)
  
  data <- data[, slotName]
  names(data) <- cellNames
  }
  
  df <- as.data.frame(table(data))
  
  # namethe first variable. 
  names(df)[1] <- slotName
  rownames(df) <- df[,1]
  
  # reorder: 
  if (decreasing) {df <- df[order(df$Freq, decreasing = TRUE),]}
  if (top) {df <- df[1:top, ]}
  
  if (return.cells) {
    cellname.list <- lapply(df[,1], function(x) cellNames[which(data == x)] ) 
    names(cellname.list) <- df[,1]
      return(cellname.list)
  }
  
  return(df)
}