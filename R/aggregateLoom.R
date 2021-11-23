# a script function to agragate several loom files into one. 
# the input should be 
# 1. read.loom object(s) (ldat files after "read.loom.matrices")
# 2. the prefix to add as cell orgin ID (or dataset ID)
# output should a loom file with 3 vectors. With cell names transformed to "Ori.ID_BC". 

aggregateLoom <- function(loomObj, Ori.ID){
  # first, check if 1. all the loomObj are loom data files; 2. the length of loom corresponds to length of Ori.ID. 
  if (! (unique((sapply(loomObj, function(x) class(x)))) == "dgCMatrix")) stop("At least one of your file is not loom data file (not class dgCMatrix). [step 1]")
  if (! (all(unique(names(sapply(loomObj, function(x) class(x)))) == c("spliced", "unspliced","ambiguous" ))))
            stop(stop("At least one of your file is not loom data file (not containing spliced, unspliced, ambiguous). [step 1]"))
  if ( !(length(loomObj) == 3*length(Ori.ID))) stop("The number of prefixs you gave is different from the number of read.loom objects. [step 1]")
  
  n <- length(Ori.ID)
  
  # fill a 2-cell loom to give formate. # here I cannot find a way to create a empty loom file.  
  loom.spliced <- loomObj[[1]][,1:2]
  loom.unspliced <- loomObj[[2]][,1:2]
  loom.ambiguous <- loomObj[[3]][,1:2]

  for(i in 0:(n-1)) 
    {
  
    # loom.obj.spliced <- list()
    # loom.obj.unspliced <- list()
    # loom.obj.ambiguous <- list()

    # below the codes lean the cell names: for example: 
    # from ‘20170630_CellRanger_F247_Dim_ScRNA_03_2:AAGGCAGTCCACGAATx’ to ‘AAGGCAGTCCACGAATx’. 
            loom.obj.spliced <- lapply(loomObj[i*3+1],function(x) {
              colnames(x) <-  gsub("_unique.bam","",gsub(".*:","",colnames(x)))
              x
            })
            loom.obj.unspliced <- lapply(loomObj[i*3+2],function(x) {
              colnames(x) <-  gsub("_unique.bam","",gsub(".*:","",colnames(x)))
              x
            })
            loom.obj.ambiguous <- lapply(loomObj[i*3+3],function(x) {
              colnames(x) <-  gsub("_unique.bam","",gsub(".*:","",colnames(x)))
              x
            })
    
    # following codes add sample prefix to every cell name, and remove "x" in the end. 
    # from ‘AAGGCAGTCCACGAATx’ to ‘ssIM1_AAGGCAGTCCACGAAT’    
            colnames(loom.obj.spliced[[1]]) <- as.vector(sapply(colnames(loom.obj.spliced[[1]]), 
                function(x) 
                  if (Ori.ID[i+1] == "") {paste(c(Ori.ID[i+1], strsplit(x, split = "x")[1]), collapse = "")} 
                else {paste(c(Ori.ID[i+1], strsplit(x, split = "x")[1]), collapse = "_")}
                ))
            
            colnames(loom.obj.unspliced[[1]]) <- as.vector(sapply(colnames(loom.obj.unspliced[[1]]), function(x) 
                  if (Ori.ID[i+1] == "") {paste(c(Ori.ID[i+1], strsplit(x, split = "x")[1]), collapse = "") } 
                  else {paste(c(Ori.ID[i+1], strsplit(x, split = "x")[1]), collapse = "_")}
                ))
            colnames(loom.obj.ambiguous[[1]]) <- as.vector(sapply(colnames(loom.obj.ambiguous[[1]]), function(x) 
              if (Ori.ID[i+1] == "") {paste(c(Ori.ID[i+1], strsplit(x, split = "x")[1]), collapse = "")} 
              else {paste(c(Ori.ID[i+1], strsplit(x, split = "x")[1]), collapse = "_")}))

    # fill lists: 
            loom.spliced <- cbind(loom.spliced,loom.obj.spliced[[1]])
            loom.unspliced <- cbind(loom.unspliced,loom.obj.unspliced[[1]])
            loom.ambiguous <- cbind(loom.ambiguous, loom.obj.ambiguous[[1]])
  }
  
  # remove the 2 cells add at beginnig of creattion. 
  loom.spliced <- loom.spliced[,-(1:2)]
  loom.unspliced <- loom.unspliced[,-(1:2)]
  loom.ambiguous <- loom.ambiguous[,-(1:2)]
  
  # combine to a loom file :
  loom <- list()
  loom <- list(spliced = loom.spliced, unspliced = loom.unspliced, ambiguous = loom.ambiguous)
  return(loom)
}