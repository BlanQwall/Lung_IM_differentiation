## filterLoom should subset a loom file with selected gene list or cell list. 


# notice: 
# 1. this function is expected to treat only one loom file at once. 
# 2. geneList and CellList should be char vectors. 

filterLoom <- function(loomObj, geneList=NULL, cellList=NULL) {
  
  # first, check if 1. all the loomObj are loom data files; 2. the length of loom corresponds to length of Ori.ID. 
  if (! (unique((sapply(loomObj, function(x) class(x)))) == "dgCMatrix")) stop("At least one of your file is not loom data file (not class dgCMatrix). [step 1]")
  if (! (all(unique(names(sapply(loomObj, function(x) class(x)))) == c("spliced", "unspliced","ambiguous" ))))
    stop(stop("At least one of your file is not loom data file (not containing spliced, unspliced, ambiguous). [step 1]"))
  
  if(is.null(geneList)) geneList <- rownames(loomObj[[1]])
  if(is.null(cellList)) cellList <- colnames(loomObj[[1]])
  
  geneList <- as.character(geneList)
  cellList <- as.character(cellList)
  
  # before all start, check all the list in seurat is also in loom: 
  if (!(all(geneList %in% rownames(loomObj[[1]])))) {
    print("Following genes are not in the loom gene list:")
    notInLoom.genes <- geneList[!(geneList %in% rownames(loomObj[[1]]))]
    print(notInLoom.genes)
    warning("No filter has been carried on because some genes cannot be found.")
    return(notInLoom.genes)
  }
  
  # and check all the cells are in loom: 
  if (!(all(cellList %in% colnames(loomObj[[1]])))) {
    print("Following cells are not in the loom cell list:")
    notInLoom.cells <- cellList[!(cellList %in% colnames(loomObj[[1]]))]
    print(notInLoom.cells)
    warning("No filter has been carried on because some cells cannot be found.")
    return(notInLoom.cells)
  }
  
  loom.obj.spliced <- loomObj[["spliced"]]
  loom.obj.unspliced <- loomObj[["unspliced"]]
  loom.obj.ambiguous <- loomObj[["ambiguous"]]
  
  loom.spliced <- loom.obj.spliced[ geneList , cellList ]
  loom.unspliced <- loom.obj.unspliced[ geneList , cellList ]
  loom.ambiguous <- loom.obj.ambiguous[ geneList , cellList ]
  
  loom <- list (spliced = loom.spliced, 
                unspliced = loom.unspliced,
                ambiguous = loom.ambiguous)
  
  return(loom)
}