seurat.setup <- function(
  path.10x, 
  merge.objects = FALSE, # if several paths given, one can choose to merge before qc. 
  prefix = NULL, # prefix before cellnames when merge. 
  project = "SeuratQC", 
  mt.percentage = 10, # percentage of mt in mt control.
  dimensionality = 1:20, 
  FindClusters.res = 0.5,
  human=FALSE
)
{
  require(Seurat)
  require(dplyr)
  if (is.null(path.10x)) {stop("path.10x not defined: a path to filtered_feature_bc_matrix/")}
  
  
  if (isFALSE(merge.objects) & length(path.10x) == 1) {
    seuratObject <- CreateSeuratObject(counts = Read10X(path.10x), 
                                     project = project, 
                                     min.cells = 3, 
                                     min.features = 200)
  
  } else {
    if (is.null(prefix) & !(length(prefix) == length(path.10x))) {stop("prefix length not equal to path length: should be 1 to 1.")} 
    object.list <- lapply(path.10x, function(x) {
      temp.object <- CreateSeuratObject(counts = Read10X(x), 
                         project = "TumorVac",
                         min.cells = 3, 
                         min.features = 200)
      index.order <- match(x, path.10x)
      levels(temp.object@meta.data$orig.ident) <- prefix[index.order]
      return(temp.object)
    })

    seuratObject <- merge(x = object.list[[1]], 
                          y = object.list[2:length(object.list)], 
                          add.cell.ids = prefix, 
                          project = project)
  }
  
  # now start qc: 
  if(human){
    seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")
  } else {
    seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^mt-")}
  
  
  # Violin plots: p1 
  p1 <- VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # feature scatter: p2
  p2 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  expr <- FetchData(object = seuratObject, vars = "percent.mt")
  seuratObject <- seuratObject[, which(x = expr < mt.percentage)]
  # seuratObject <- subset(seuratObject, cells = WhichCells(seuratObject, expression = 'percent.mt' < mt.percentage)) # mt.pct change to 10
  
  # normalization: 
  seuratObject <- NormalizeData(seuratObject)
  
  ## Identification of highly variable features (feature selection)
  seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top20 <- head(VariableFeatures(seuratObject), 20)
  
  # plot variable features with labels: p3
  p3 <- LabelPoints(plot = VariableFeaturePlot(seuratObject), points = top20, repel = TRUE)
  
  ##Scaling the data
  seuratObject <- ScaleData(seuratObject, features = rownames(seuratObject))
  
  # PCA: 
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  
  # plot PCA: p4
  p4 <- DimPlot(seuratObject, reduction = "pca", group.by = "orig.ident")
  
  ## Determine the ‘dimensionality’ of the dataset
  seuratObject <- JackStraw(seuratObject, num.replicate = 100)
  seuratObject <- ScoreJackStraw(seuratObject, dims = dimensionality)
  
  # jackstrawPlot: p5
  p5 <- JackStrawPlot(seuratObject, dims = dimensionality)
  
  # ElbowPlot: p6
  p6 <- ElbowPlot(seuratObject)
  
  
  seuratObject <- FindNeighbors(seuratObject, dims = dimensionality)
  seuratObject <- FindClusters(seuratObject, resolution = FindClusters.res)
  
  ## Run non-linear dimensional reduction (UMAP/tSNE)
  seuratObject <- RunUMAP(seuratObject, dims = dimensionality)
  
  # UMAP plot: p7
  p7 <- DimPlot(seuratObject, reduction = "umap")
  
  seuratObject <- RunTSNE(seuratObject, dims = dimensionality)
  
  # TSNE plot: p8
  p8 <- DimPlot(seuratObject, reduction = "tsne")
  
  result.list <- list(seuratObject = seuratObject,
                      plots = list(feature_vln = p1, 
                                   RNA_mt.pct_scatter = p2, 
                                   variable_features = p3, 
                                   PCA_plot = p4,
                                   JackStrawPlot = p5, 
                                   ElbowPlot = p6, 
                                   UMAP_plot = p7, 
                                   TSNE_plot = p8))
  
  return(result.list)
  
}