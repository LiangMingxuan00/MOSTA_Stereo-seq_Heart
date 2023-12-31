
my_convert <- function(h5adfile,rdsfile,raw_coord_path){
	infile <- h5adfile
	outfile <- rdsfile
	Convert(infile, dest = "h5seurat", assay = "Spatial", overwrite = TRUE)
	# convert h5ad as h5seurat, which means a seurat-object format stored in h5

	h5file <- paste(paste(unlist(strsplit(infile, "h5ad", fixed = TRUE)), collapse='h5ad'), "h5seurat", sep="")
	print(paste(c("Finished! Converting h5ad to h5seurat file at:", h5file), sep=" ", collapse=NULL))

	object <- LoadH5Seurat(h5file)
	print(paste(c("Successfully load h5seurat:", h5file), sep=" ", collapse=NULL))

	# spatial already transform to `Spatial`` in assays
	if (!is.null(object@reductions$spatial)) {
	  object@reductions$spatial <- NULL
	}

	# convert stereopy SCT result to seurat SCT result
	if (
	    !is.null(object@misc$sct_counts) &&
	    !is.null(object@misc$sct_data) &&
	    !is.null(object@misc$sct_scale) &&
	    !is.null(object@misc$sct_cellname) &&
	    !is.null(object@misc$sct_genename) &&
	    !is.null(object@misc$sct_top_features)
	  ) {
	  sct.assay.out <- CreateAssayObject(counts=object[['Spatial']]@counts, check.matrix=FALSE)
	  sct.assay.out <- SetAssayData(
	      object = sct.assay.out,
	      slot = "data",
	      new.data = log1p(x=GetAssayData(object=sct.assay.out, slot="counts"))
	    )
	  sct.assay.out@scale.data <- as.matrix(object@misc$sct_scale)
	  colnames(sct.assay.out@scale.data) <- object@misc$sct_cellname
	  rownames(sct.assay.out@scale.data) <- object@misc$sct_scale_genename
	  sct.assay.out <- Seurat:::SCTAssay(sct.assay.out, assay.orig='Spatial')
	  Seurat::VariableFeatures(object = sct.assay.out) <- object@misc$sct_top_features
	  object[['SCT']] <- sct.assay.out
	  DefaultAssay(object=object) <- 'SCT'

	  # TODO: tag the reductions as SCT, this will influence the find_cluster choice of data
	  object@reductions$pca@assay.used <- 'SCT'
	  object@reductions$umap@assay.used <- 'SCT'
	  assay.used <- 'SCT'
	  print("Finished! Got SCTransform result in object, create a new SCTAssay and set it as default assay.")
	} else {
	  # TODO: we now only save raw counts, try not to add raw counts to .data, do NormalizeData to fit this
	  object <- NormalizeData(object)
	  assay.used <- 'Spatial'
	  print("Finished! Got raw counts only, auto create log-normalize data.")
	}


	# TODO follow with old code, don't touch
	print("Start add image...This may take some minutes...(~.~)")
	# add image
	raw_coord <- read.csv(raw_coord_path)
	colnames(raw_coord) <- c("x","y")
	cell_coords <- unique(raw_coord[, c('x', 'y')])
	cell_coords['cell'] <- row.names(object@meta.data)
	cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
	cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1

	# object of images$slice1@image, all illustrated as 1 since no concrete pic
	tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))

	# object of images$slice1@coordinates, concrete coordinate of X and Y
	tissue_positions_list <- data.frame(row.names = cell_coords$cell,
	                                    tissue = 1,
	                                    row = cell_coords$y, col = cell_coords$x,
	                                    imagerow = cell_coords$y, imagecol = cell_coords$x)
	# @images$slice1@scale.factors
	scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1, tissue_hires_scalef = 1, tissue_lowres_scalef = 1))

	# generate object @images$slice1
	generate_BGI_spatial <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
	  if (filter.matrix) {
	    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
	  }
	  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
	  spot.radius <- unnormalized.radius / max(dim(x = image))
	  return(new(Class = 'VisiumV1',
	             image = image,
	             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
	                                          fiducial = scale.factors$fiducial_diameter_fullres,
	                                          hires = scale.factors$tissue_hires_scalef,
	                                          lowres = scale.factors$tissue_lowres_scalef),
	             coordinates = tissue.positions,
	             spot.radius = spot.radius))
	}

	BGI_spatial <- generate_BGI_spatial(image = tissue_lowres_image,
	                                    scale.factors = fromJSON(scalefactors_json),
	                                    tissue.positions = tissue_positions_list)

	# can be thought of as a background of spatial
	# import image into seurat object
	object@images[['slice1']] <- BGI_spatial
	object@images$slice1@key <- "slice1_"
	object@images$slice1@assay <- assay.used


	print("Finished add image...Start to saveRDS...")
	saveRDS(object, outfile)

}