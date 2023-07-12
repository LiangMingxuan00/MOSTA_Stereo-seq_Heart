library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

store.path <- c("//biomedja01/disk1/lmx_home/thesis/processed_data/")
plot_name <- c("")

saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
	ggsave(filename, plot, path = store.path, width = width, height = height)
}



import anndata
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame

adata = anndata.read_h5ad("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.h5ad")

adata.obsm['spatial']

bdata = adata[adata.obs.cell_type == "B"]
bdata

bdata = adata[adata.obs.timepoint == "E9.5"]

spatialinfo  = bdata.obsm['spatial']

c = np.array(spatialinfo)

d = DataFrame(spatialinfo)

d.to_csv("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/e9_5spatial.csv",sep=',')




adata = ReadH5AD("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.h5ad")








library(Seurat)
library(SeuratData)
library(SeuratDisk)
# Converting from AnnData to Seurat via h5Seurat

Convert("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.h5ad", dest="h5seurat", overwrite=T)
mosta <- LoadH5Seurat("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.h5seurat")

saveRDS(mosta,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.rds")

mosta <- readRDS("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.rds")

mosta.list <- SplitObject(mosta, split.by = "timepoint")



saveRDS(mosta.list$E9.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E9.5.rds")
saveRDS(mosta.list$E10.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E10.5.rds")
saveRDS(mosta.list$E11.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E11.5.rds")
saveRDS(mosta.list$E12.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E12.5.rds")
saveRDS(mosta.list$E13.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E13.5.rds")
saveRDS(mosta.list$E14.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E14.5.rds")
saveRDS(mosta.list$E15.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E15.5.rds")
saveRDS(mosta.list$E16.5,"//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E16.5.rds")





a <- readRDS("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_E9.5.rds")



data <- read.csv("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/e9_5spatial.csv")

coord = data[,2:3]


barcodes <- rownames(a@meta.data)


rownames(coord) <- rownames(a@meta.data)


slide.seq = CreateSeuratObject(counts = COUNTS_MTX, assay="Spatial")

coord.df = data.frame(x=X, y=Y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = BARCODES





slide.seq@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )




a@images$image = new(Class = "SlideSeq",assay = "RNA", key = "image_",coordinates = coord,scale.factors=1)



a@images$image <- new(
  Class = 'VisiumV1',
  assay = "RNA",
  key = "image_",
  scale.factors = 1, # scale factors
  coordinates = coord, # data frame with rows as spots and columns as axes
)




