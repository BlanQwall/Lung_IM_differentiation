# this function will return a list containing: $integrated.seuratObject and $plots

seurat.make.integrated <- function(
  seuratObjects, # input: c(object1, object2...)
  prefix = NULL, # this will give a name to metadata named "object_before_integrated". 
  # findVariableFeatures.nfeature = 2000, # cannot use. 
  dimensionality = 1:30, 
  SCTransf = TRUE, # SCTransformation will be used if true (see Seurat SCTransform)
  ...
)
{
  require(Seurat)
  
  if (is.null(seuratObjects) | length(seuratObjects) < 2) {stop("seuratObjects must be defined as multiple.")}
  if (is.null(prefix)) {stop("prefix not given. ")}
  if (!(length(prefix) == length(seuratObjects) )) { stop("prefix length doesn't equal to seuratObjects length")}
  
  list.so <- as.list(seuratObjects)
  names(list.so) <- prefix
  
  
  # give prefix name to metadata named "object_before_integrated". 
  for (i in 1:length(list.so)) {
    list.so[[i]][["object_before_integrated"]] <- prefix[i]
  }
  
  # make the RNA as the default assay: in case of integrated analysis oreadly done, this will issue error if not not using "RNA. 
  for (i in 1:length(list.so)) {
    DefaultAssay(list.so[[i]]) <- "RNA"
  }
  
  # return(list.so)
  # 
  
  # standard process: preprocessing (log-normalization), and identify variable features individually for each.
  list.so <- lapply(list.so, function(x) {
    if (SCTransf) {
      suppressWarnings(
                x <- SCTransform(x, verbose = FALSE)
      )
    } else {    x <- NormalizeData(x, verbose = FALSE)
                x <- FindVariableFeatures(x, selection.method = "vst",
                              nfeatures = 2000, # FindVariableFeatures.nfeature is not usable. 
                              verbose = FALSE)
      }
    return(x)
  } )

  # integration:
  if (SCTransf) {
    integrateion.features <- SelectIntegrationFeatures(object.list = list.so, nfeatures = 3000)
    list.so <- PrepSCTIntegration(object.list = list.so, anchor.features = integrateion.features, 
                                        verbose = FALSE)
    anchors <- FindIntegrationAnchors(object.list = list.so, normalization.method = "SCT", 
                                               anchor.features = integrateion.features, verbose = FALSE, ...)
    integrated.so <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                                         verbose = FALSE)
  } else {
    anchors <- FindIntegrationAnchors(object.list = list.so, dims = dimensionality, ...)
    integrated.so <- IntegrateData(anchorset = anchors)
  }



  # cluster
  require(ggplot2)

  # set default assay:
  DefaultAssay(integrated.so) <- "integrated"

  # Run the standard workflow for visualization and clustering
  integrated.so <- ScaleData(integrated.so, verbose = FALSE)
  integrated.so <- RunPCA(integrated.so,
                          npcs = tail(dimensionality, n = 1),  # use the last number of dimensionality to make pca.
                          verbose = FALSE)
  integrated.so <- RunUMAP(integrated.so, reduction = "pca", dims = dimensionality)
  integrated.so <- RunTSNE(integrated.so, reduction = "pca", dims = dimensionality)
  integrated.so <- FindNeighbors(integrated.so, reduction = "pca", dims = dimensionality)

  # plot pca

  p.pca <- DimPlot(integrated.so, reduction = "pca")
  p.pca.group <- DimPlot(integrated.so, reduction = "pca",
                         split.by = "object_before_integrated",
                         label = TRUE,
                         repel = TRUE) + NoLegend()

  # TSNE plot:

  p.tsne <- DimPlot(integrated.so, reduction = "tsne")
  p.tsne.group <- DimPlot(integrated.so, reduction = "tsne",
                         split.by = "object_before_integrated",
                         label = TRUE,
                         repel = TRUE) + NoLegend()

  # UMAP plot:

  p.umap <- DimPlot(integrated.so, reduction = "umap")
  p.umap.group <- DimPlot(integrated.so, reduction = "umap",
                          split.by = "object_before_integrated",
                          label = TRUE,
                          repel = TRUE) + NoLegend()


  restuls <- list(integrated.seuratObject = integrated.so,
                  plots = list (plot.pca = p.pca,
                                plot.pca.group = p.pca.group,
                                plot.tsne = p.tsne,
                                plot.tsne.group = p.tsne.group,
                                plot.umap = p.umap,
                                plot.umap.group = p.umap.group))

  return(restuls)
}
  