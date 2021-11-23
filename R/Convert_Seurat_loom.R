Convert.seurat_loom <- function( # modified from Seurat::Convert
  from,
  to,
  filename,
  overwrite = FALSE,
  display.progress = TRUE,
  anndata.raw = "raw.data",
  anndata.X = "data",
  ...
) {
  object.to <- switch(
    EXPR = to,
    'loom' = {
      if (!'loomR' %in% rownames(x = installed.packages())) {
        stop("Please install loomR from GitHub before converting to a loom object")
      }
      
      # deal with list: 
      ldat <- from$ldat
      from <- from$seurat
      cell.order.ldat <- colnames(ldat[[1]])
      gene.order.ldat <- rownames(ldat[[1]])
      
      cell.order <- colnames(from)
      gene.order <- rownames(x = GetAssayData(from, assay = "RNA"))
      
      # to filter only common genes and cells in seurat: 
      cell.order <- intersect(cell.order, cell.order.ldat)
      gene.order <- intersect(gene.order, gene.order.ldat)
      
      from <- from[gene.order, cell.order]
      ldat <- lapply(ldat, function(x) {x <- x[gene.order, cell.order]} )
      
      loomfile <- loomR::create(
        filename = filename,
        data = GetAssayData(from, assay = "RNA", slot = "counts")[gene.order, cell.order],
        cell.attrs = from@meta.data[cell.order, ],
        # layers = list('norm_data' = t(x = as.matrix(GetAssayData(from, assay = "RNA")[, cell.order]))),
        layers = lapply(ldat, function(x) t(as.matrix(x))),
        chunk.size = NULL,
        chunk.dims = NULL,
        overwrite = overwrite,
        display.progress = display.progress
      )
      if (nrow(x = from@assays$RNA@meta.features) > 0 & ncol(x = from@assays$RNA@meta.features) > 0) {
        hvg.info <- from@assays$RNA@meta.features
        colnames(x = hvg.info) <- gsub(
          pattern = '.',
          replacement = '_',
          x = colnames(x = hvg.info),
          fixed = TRUE
        )
        loomfile$add.row.attribute(hvg.info[gene.order, ])
      }
      if (length(x = from@assays$RNA@var.features) > 0) {
        loomfile$add.row.attribute(list('var_genes' = gene.order %in% from@assays$RNA@var.features))
      }
      # if (!is.null(x = from@assays$RNA@scale.data) && dim(x = from@assays$RNA@scale.data) != c(1, 1)) {
      #   loomfile$add.layer(list(
      #     'scale_data' = as.matrix(x = t(x = as.data.frame(x = from@assays$RNA@scale.data)[gene.order, cell.order]))
      #   ))
      # }
      for (dim.reduc in names(x = from@reductions)) {
        cell.embeddings <- from@reductions[[dim.reduc]]@cell.embeddings
        ce.dims <- unique(x = dim(x = cell.embeddings))
        if (length(x = ce.dims) != 1 || ce.dims != 0) {
          if (nrow(x = cell.embeddings) < ncol(x = GetAssayData(from, assay = "RNA"))) {
            cell.embeddings.padded <- matrix(
              nrow = length(x = colnames(from)),
              ncol = ncol(x = cell.embeddings)
            )
            if (is.null(x = rownames(x = cell.embeddings)) || is.null(x = colnames(from))) {
              pad.order <- 1:nrow(x = cell.embeddings)
            } else {
              pad.order <- match(
                x = rownames(x = cell.embeddings),
                table = colnames(from)
              )
            }
            cell.embeddings.padded[pad.order, ] <- cell.embeddings
          } else if (nrow(x = cell.embeddings) > ncol(x = GetAssayData(from, assay = "RNA"))) {
            stop("Cannot have more cells in the dimmensional reduction than in the dataset")
          } else {
            cell.embeddings.padded <- cell.embeddings
          }
          cell.embeddings.padded <- list(cell.embeddings.padded)
          names(x = cell.embeddings.padded) <- paste0(dim.reduc, '_cell_embeddings')
          loomfile$add.col.attribute(cell.embeddings.padded)
        }
        feature.loadings <- from@reductions[[dim.reduc]]@feature.loadings.projected
        gl.dims <- unique(x = dim(x = feature.loadings))
        if (length(x = gl.dims) == 1 && gl.dims == 0) {
          feature.loadings <- from@reductions[[dim.reduc]]@feature.loadings
        }
        gl.dims <- unique(x = dim(x = feature.loadings))
        if (length(x = gl.dims) != 1 || gl.dims != 0) {
          if (nrow(x = feature.loadings) < nrow(x = GetAssayData(from, assay = "RNA"))) {
            feature.loadings.padded <- matrix(
              nrow = nrow(x = GetAssayData(from, assay = "RNA")),
              ncol = ncol(x = feature.loadings)
            )
            if (is.null(x = rownames(x = feature.loadings)) || is.null(x = rownames(x = GetAssayData(from, assay = "RNA")))) {
              pad.order <- 1:nrow(x = feature.loadings)
            } else {
              pad.order <- match(
                x = rownames(x = feature.loadings),
                table = rownames(x = GetAssayData(from, assay = "RNA"))
              )
            }
            feature.loadings.padded[pad.order, ] <- feature.loadings
          } else if (nrow(x = feature.loadings) > nrow(x = GetAssayData(from, assay = "RNA"))) {
            stop("Cannot have more genes in the dimmensional reduction than in the dataset")
          } else {
            feature.loadings.padded <- feature.loadings
          }
          feature.loadings.padded <- list(feature.loadings.padded)
          names(x = feature.loadings.padded) <- paste0(dim.reduc, '_gene_loadings')
          loomfile$add.row.attribute(feature.loadings.padded)
        }
      }
      loomfile
    },
    stop(paste0("Cannot convert Seurat objects to class '", to, "', only support convert Seurat/loom list to loom."))
  )
  return(object.to)
}
