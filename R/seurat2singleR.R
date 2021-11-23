
seurat2singleR <- function(
  seuratObject, 
  ref = "ImmGenData", 
  labels = "label.main", 
  method = "single",
  reduction = "umap"
)
  
{
  require(Seurat)
  require(SingleR)
  
  # check object version: 
  if (seuratObject@version < 3) {
    cat(paste("The current object is Seurat Object version", as.character (seuratObject@version), "\n" ))
    stop("Old version object.")
  }
  obj <- seuratObject
  # find a reduction to use: priority = umap > tsne > pca. 
  if (is.null(obj@reductions[reduction])) { reduction = "umap" }
  if (is.null(obj@reductions[reduction])) { reduction = "tsne" }
  if (is.null(obj@reductions[reduction])) { reduction = "pca" }
  if (is.null(obj@reductions[reduction])) {stop("None of reduction (UMAP, TSNE, PCA) is available.")}
  
  data <- GetAssayData(object = seuratObject, slot = "data", assay = "RNA" )
  
  # prepare for ref and lables: 
    ref = do.call(ref, args = list())
    labels = slot(ref, "colData")[, labels]
  singler = SingleR(test = data, ref = ref, labels = labels)
  
  singler$meta.data$orig.ident <- seuratObject$orig.ident
  # singler$meta.data$xy <- seuratObject@reductions[[reduction]]@cell.embeddings # use TSNE to discriminate subsets. 
  # singler$meta.data$clusters = seuratObject@active.ident # no need, we create singler object with this slot. 
  
  # cat("Saving SinglerObject.rds ... \n")
  # saveRDS(singler,file="NewSinglerObject.rds")
  # cat("\n")
  # cat("Now you can load this rds file to http://comphealth.ucsf.edu/SingleR. \n")
  return(singler)
}


# --------------- old version --------------------
# seurat2singleR <- function(
#   seuratObject, 
#   project = "Project_singleR", 
#   species = "Mouse",  # could be Human
#   reduction = "umap"
# )
# 
# {
#   require(Seurat)
#   require(SingleR)
#   # require(RColorBrewer)
#   # source("~/R/x86_64-pc-linux-gnu-library/3.6/SingleR/R/SingleR.Create.R", )
#   # source("~/R/x86_64-pc-linux-gnu-library/3.6/SingleR/R/SingleR.Object.R")
#   # source("~/R/x86_64-pc-linux-gnu-library/3.6/SingleR/R/SingleR.References.R")
#   
#   # check object version: 
#   if (seuratObject@version < 3) {
#     cat(paste("The current object is Seurat Object version", as.character (seuratObject@version), "\n" ))
#     stop("Old version object.")
#   }
#   obj <- seuratObject
#   # find a reduction to use: priority = umap > tsne > pca. 
#   if (is.null(obj@reductions[reduction])) { reduction = "umap" }
#   if (is.null(obj@reductions[reduction])) { reduction = "tsne" }
#   if (is.null(obj@reductions[reduction])) { reduction = "pca" }
#   if (is.null(obj@reductions[reduction])) {stop("None of reduction (UMAP, TSNE, PCA) is available.")}
#   
#   data <- GetAssayData(object = seuratObject, slot = "counts", assay = "RNA" )
#   singler = CreateSinglerObject(data, annot = NULL, project.name = project, min.genes = 0,
#                                 technology = "10X", species = species, citation = "",
#                                 ref.list = list(), normalize.gene.length = F, variable.genes = "de",
#                                 fine.tune = FALSE, do.signatures = T, clusters = seuratObject@active.ident, do.main.types = T, 
#                                 reduce.file.size = T, numCores = SingleR.numCores)
#   
#   singler$meta.data$orig.ident <- seuratObject$orig.ident
#   singler$meta.data$xy <- seuratObject@reductions[[reduction]]@cell.embeddings # use TSNE to discriminate subsets. 
#   singler$meta.data$clusters = seuratObject@active.ident # no need, we create singler object with this slot. 
#   
#   # see update from https://github.com/dviraran/SingleR
#   
#   singler.new = convertSingleR2Browser(singler)
#   cat("This function will also create a NewSinglerObject.rds file for web use. \n")
#   cat("Saving NewSinglerObject.rds ... \n")
#   saveRDS(singler.new,file="NewSinglerObject.rds")
#   cat("\n")
#   cat("Now you can load this rds file to http://comphealth.ucsf.edu/SingleR. \n")
# 
# }