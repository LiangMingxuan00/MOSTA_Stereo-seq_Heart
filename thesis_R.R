### Thesis: Integration scRNA-seq and ST for embroyonic heart development
### Author: Liang Mingxuan
### Student ID: 3036082432


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Python script
import anndata
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns


# Requirements should be satisfied by a PEP 517 installer. If you are using pip, you can try `pip install --use-pep517`.
# Using: "pip install --use-pep517 stereopy"
import stereo as st


mouse_all_stage = anndata.read_h5ad("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/Mouse_embryo_all_stage.h5ad") #load all stage mouse embroyo data

mouse_E9_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E9.5"]
mouse_E10_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E10.5"]
mouse_E11_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E11.5"]
mouse_E12_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E12.5"]
mouse_E13_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E13.5"]
mouse_E14_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E14.5"]
mouse_E15_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E15.5"]
mouse_E16_5 = mouse_all_stage[mouse_all_stage.obs.timepoint == "E16.5"]

mouse_E9_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E9_5.h5ad")
mouse_E10_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E10_5.h5ad")
mouse_E11_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E11_5.h5ad")
mouse_E12_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E12_5.h5ad")
mouse_E13_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E13_5.h5ad")
mouse_E14_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E14_5.h5ad")
mouse_E15_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E15_5.h5ad")
mouse_E16_5.write("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E16_5.h5ad")




###output the spatial information
### spatialinfo(): extract the spatial localization for the convert h5ad to seurat
import numpy as np
import pandas as pd

def spatialinfo(bdata, filepath):
    spatial_data = np.array(bdata.obsm['spatial'])
    df = pd.DataFrame(spatial_data)
    df.to_csv(filepath, sep=',', index=False)

spatialinfo(mouse_E9_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E9_5spatial.csv")
spatialinfo(mouse_E10_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E10_5spatial.csv")
spatialinfo(mouse_E11_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E11_5spatial.csv")
spatialinfo(mouse_E12_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E12_5spatial.csv")
spatialinfo(mouse_E13_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E13_5spatial.csv")
spatialinfo(mouse_E14_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E14_5spatial.csv")
spatialinfo(mouse_E15_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E15_5spatial.csv")
spatialinfo(mouse_E16_5, "//biomedja01/disk1/lmx_home/thesis/bgi_dataset/E16_5spatial.csv")



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### R script
### Convert h5ad to seurat object
### SeuratDisk local enviorment "cellchat"
library(dplyr)
library(rjson)
library(Seurat)
library(ggplot2)
library(SeuratDisk)
source("/biomedja01/disk1/lmx_home/thesis/thesis_myconvert.R") ### import the convert function


namelist <- c("E9_5","E10_5","E11_5","E12_5","E13_5","E14_5","E15_5","E16_5")

for(i in namelist){
  inputfile <- c(paste0("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_",i,".h5ad"))
  outputfile <- c(paste0("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_",i,".rds"))
  csvfile <- c(paste0("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/",i,"spatial.csv"))
  print(inputfile)
  print(outputfile)
  print(csvfile)
  my_convert(inputfile,outputfile,csvfile) ### convert function
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### R script
### Basic analysis of dataset
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggalluvial)
library(tidydr)
library(paletteer)

###saveplot(): save the pdf file
saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
  ggsave(filename, plot, path = store.path, width = width, height = height,dpi=600)
}
store.path <- c("/biomedja01/disk1/lmx_home/thesis/processed_data")
plot_name <- c("")

options(future.globals.maxSize = 100000 * 1024^2) ###Set CPU


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#




### loop for the list to load the data in seurat
namelist <- c("E9_5","E10_5","E11_5","E12_5","E13_5","E14_5","E15_5","E16_5")
for(i in namelist){
  assign(i, readRDS(paste0("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_", i, ".rds")))
}


##
E13_5 <- readRDS("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E13_5.rds")
E14_5 <- readRDS("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E14_5.rds")
E15_5 <- readRDS("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E15_5.rds")
E16_5 <- readRDS("//biomedja01/disk1/lmx_home/thesis/bgi_dataset/mouse_E16_5.rds")


### draw the orginal plot of each stage sample
varlist <- list(E9_5= E9_5,E10_5=E10_5,E10_5=E10_5,E11_5=E11_5,E12_5=E12_5,E13_5=E13_5,E14_5=E14_5,E15_5=E15_5,E16_5=E16_5)
x=1

for(i in varlist){
  SpatialSample <- SpatialDimPlot(i, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
  name <- paste0("Original_Spatialplot_",names(varlist[x]))
  saveplot(plot_name,SpatialSample,name,store.path,24,24)
  x= x+1
}



### standard workflow process the object
reduction_sample <- function(sample) {
  sample  <- SCTransform(sample, assay = "Spatial", verbose = FALSE)
  sample  <- RunPCA(sample , assay = "SCT", verbose = F)
  sample  <- FindNeighbors(sample , reduction = "pca", dims = 1:30, verbose = F)
  sample  <- FindClusters(sample , resolution = 0.5 ,verbose = F,)
  sample  <- RunUMAP(sample , dims = 1:30, verbose = F)
  return(sample)
}



E9_5 <- reduction_sample(E9_5)
E10_5 <- reduction_sample(E10_5)
E11_5 <- reduction_sample(E11_5)
E12_5 <- reduction_sample(E12_5)
E13_5 <- reduction_sample(E13_5)
E14_5 <- reduction_sample(E14_5)
E15_5 <- reduction_sample(E15_5)
E16_5 <- reduction_sample(E16_5)

saveRDS(E9_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5.rds")
saveRDS(E10_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5.rds")
saveRDS(E11_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5.rds")
saveRDS(E12_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5.rds")
saveRDS(E13_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E13_5.rds")
saveRDS(E14_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E14_5.rds")
saveRDS(E15_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E15_5.rds")
saveRDS(E16_5,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E16_5.rds")


E9_5  <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5.rds")
E10_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5.rds")
E11_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5.rds")
E12_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5.rds")
E13_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E13_5.rds")
E14_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E14_5.rds")
E15_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E15_5.rds")
E16_5 <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E16_5.rds")


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
###Using Heart Marker to subset the Heart cluster
Dotplot_marker <- function(object,markerlist,name,width,height){
  
  genes_to_check = list(Heart = markerlist)

  p_all_markers = DotPlot(object, 
                      features = genes_to_check,
                      scale = T,assay='SCT',group.by = "SCT_snn_res.0.5" )+
    theme_bw()+
    scale_color_continuous(low="grey",high =  "red")+
    theme(legend.position = "right",legend.box = "vertical",
    legend.margin=margin(t= 0, unit='cm'),
    #legend.spacing = unit(0,"in"),
    axis.text.x  = element_text(color="black",size=12,angle = 45,vjust = 0.5, hjust=0.5),
    axis.text.y  = element_text(color="black",size=12),
    legend.text = element_text(size =12,color="black"),
    legend.title = element_text(size =12,color="black")
    ) 

    saveplot(plot_name,p_all_markers,name,store.path,width,height)
}



markerlist <- c("Actc1","Tnnt2","Myl4","Myl7","Tnni1","Myh6")
Dotplot_marker(E9_5,markerlist,"Heart_Dotplot_E9_5",8,4)
Dotplot_marker(E10_5,markerlist,"Heart_Dotplot_E10_5",8,4)
Dotplot_marker(E11_5,markerlist,"Heart_Dotplot_E11_5",8,6)
Dotplot_marker(E12_5,markerlist,"Heart_Dotplot_E12_5",8,7)
Dotplot_marker(E13_5,markerlist,"Heart_Dotplot_E13_5",8,7)
Dotplot_marker(E14_5,markerlist,"Heart_Dotplot_E14_5",8,7)
Dotplot_marker(E15_5,markerlist,"Heart_Dotplot_E15_5",8,7)
Dotplot_marker(E16_5,markerlist,"Heart_Dotplot_E16_5",8,7)




### observe and subset Heart cluster based on Dotplot

### E9.5
E9_5@meta.data$new_ident <- E9_5@meta.data$SCT_snn_res.0.5
E9_5@meta.data$new_ident <- "Other"
E9_5@meta.data[E9_5@meta.data$SCT_snn_res.0.5=="5",]$new_ident <- "Heart"
Idents(E9_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E9_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E9_5_Heart",store.path,24,24)

### E10.5
E10_5@meta.data$new_ident <- E10_5@meta.data$SCT_snn_res.0.5
E10_5@meta.data$new_ident <- "Other"
E10_5@meta.data[E10_5@meta.data$SCT_snn_res.0.5=="5",]$new_ident <- "Heart"
Idents(E10_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E10_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E10_5_Heart",store.path,24,24)

### E11.5
E11_5@meta.data$new_ident <- E11_5@meta.data$SCT_snn_res.0.5
E11_5@meta.data$new_ident <- "Other"
E11_5@meta.data[E11_5@meta.data$SCT_snn_res.0.5=="5",]$new_ident <- "Heart"
E11_5@meta.data[E11_5@meta.data$SCT_snn_res.0.5=="22",]$new_ident <- "Heart"
Idents(E11_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E11_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E11_5_Heart",store.path,24,24)


### E12.5
E12_5@meta.data$new_ident <- E12_5@meta.data$SCT_snn_res.0.5
E12_5@meta.data$new_ident <- "Other"
E12_5@meta.data[E12_5@meta.data$SCT_snn_res.0.5=="9",]$new_ident <- "Heart"
Idents(E12_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E12_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E12_5_Heart",store.path,48,48)

### E13.5
E13_5@meta.data$new_ident <- E13_5@meta.data$SCT_snn_res.0.5
E13_5@meta.data$new_ident <- "Other"
# E13_5@meta.data[E13_5@meta.data$SCT_snn_res.0.5=="9",]$new_ident <- "Heart"
E13_5@meta.data[E13_5@meta.data$SCT_snn_res.0.5=="8",]$new_ident <- "Heart"
Idents(E13_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E13_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E13_5_Heart",store.path,24,24)

### E14.5
E14_5@meta.data$new_ident <- E14_5@meta.data$SCT_snn_res.0.5
E14_5@meta.data$new_ident <- "Other"
# E14_5@meta.data[E14_5@meta.data$SCT_snn_res.0.5=="0",]$new_ident <- "Heart"
E14_5@meta.data[E14_5@meta.data$SCT_snn_res.0.5=="12",]$new_ident <- "Heart"
Idents(E14_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E14_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E14_5_Heart",store.path,24,24)

### E15.5
E15_5@meta.data$new_ident <- E15_5@meta.data$SCT_snn_res.0.5
E15_5@meta.data$new_ident <- "Other"
# E15_5@meta.data[E15_5@meta.data$SCT_snn_res.0.5=="3",]$new_ident <- "Heart"
# E15_5@meta.data[E15_5@meta.data$SCT_snn_res.0.5=="28",]$new_ident <- "Heart"
E15_5@meta.data[E15_5@meta.data$SCT_snn_res.0.5=="10",]$new_ident <- "Heart"
Idents(E15_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E15_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E15_5_Heart",store.path,24,24)

### E16.5
E16_5@meta.data$new_ident <- E16_5@meta.data$SCT_snn_res.0.5
E16_5@meta.data$new_ident <- "Other"
E16_5@meta.data[E16_5@meta.data$SCT_snn_res.0.5=="24",]$new_ident <- "Heart"
E16_5@meta.data[E16_5@meta.data$SCT_snn_res.0.5=="9",]$new_ident <- "Heart"
Idents(E16_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E16_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E16_5_Heart",store.path,24,24)



### subset the heart cluster
Idents(E9_5) <- "SCT_snn_res.0.5"
Idents(E10_5) <- "SCT_snn_res.0.5"
Idents(E11_5) <- "SCT_snn_res.0.5"
Idents(E12_5) <- "SCT_snn_res.0.5"
Idents(E13_5) <- "SCT_snn_res.0.5"
Idents(E14_5) <- "SCT_snn_res.0.5"
Idents(E15_5) <- "SCT_snn_res.0.5"
Idents(E16_5) <- "SCT_snn_res.0.5"

E9_5_heart <- subset(E9_5,idents=c("5"))
E10_5_heart <- subset(E10_5,idents=c("5"))
E11_5_heart <- subset(E11_5,idents=c("5","22"))
E12_5_heart <- subset(E12_5,idents=c("9"))
E13_5_heart <- subset(E13_5,idents=c("8"))
E14_5_heart <- subset(E14_5,idents=c("12"))
E15_5_heart <- subset(E15_5,idents=c("10"))
E16_5_heart <- subset(E16_5,idents=c("9","24"))


############################################################################################################################################
### test code using annotation
# Idents(E9_5) <- "annotation"
# Idents(E10_5) <- "annotation"
# Idents(E11_5) <- "annotation"
# Idents(E12_5) <- "annotation"
# Idents(E13_5) <- "annotation"
# Idents(E14_5) <- "annotation"
# Idents(E15_5) <- "annotation"
# Idents(E16_5) <- "annotation"

# E9_5_heart <- subset(E9_5,idents = c("Heart"))
# E10_5_heart <- subset(E10_5,idents = c("Heart"))
# E11_5_heart <- subset(E11_5,idents = c("Heart"))
# E12_5_heart <- subset(E12_5,idents = c("Heart"))
# E13_5_heart <- subset(E13_5,idents = c("Heart"))
# E14_5_heart <- subset(E14_5,idents = c("Heart"))
# E15_5_heart <- subset(E15_5,idents = c("Heart"))
# E16_5_heart <- subset(E16_5,idents = c("Heart"))


# SpatialSample <- SpatialDimPlot(E9_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E9_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E10_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E10_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E11_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E11_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E12_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E12_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E13_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E13_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E14_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E14_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E15_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E15_5",store.path,24,24)

# SpatialSample <- SpatialDimPlot(E16_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
# saveplot(plot_name,SpatialSample,"Spatialplot_Heart_test_E16_5",store.path,24,24)
##########################################################################################################################################################################



### Based on histolgy to remove some cell away from the main part of heart
### Based on image

E11_5_heart_subset <- subset(E11_5_heart, slice1_imagecol>160,invert = TRUE)
SpatialSample <- SpatialDimPlot(E11_5_heart_subset, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_E11_5_Heart_subset",store.path,24,24)


saveRDS(E9_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5_heart.rds")
saveRDS(E10_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5_heart.rds")
saveRDS(E11_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5_heart.rds")
saveRDS(E11_5_heart_subset,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5_heart_subset.rds")
saveRDS(E12_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5_heart.rds")
saveRDS(E13_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E13_5_heart.rds")
saveRDS(E14_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E14_5_heart.rds")
saveRDS(E15_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E15_5_heart.rds")
saveRDS(E16_5_heart,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E16_5_heart.rds")




#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
###R script
### In this part, we will integration all stage heart sample together to identify the cell types

E9_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5_heart.rds")
E10_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5_heart.rds")
E11_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5_heart_subset.rds") ### subsetversion
E12_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5_heart.rds")
E13_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E13_5_heart.rds")
E14_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E14_5_heart.rds")
E15_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E15_5_heart.rds")
E16_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E16_5_heart.rds")


SpatialSample <- SpatialDimPlot(E9_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E9_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E10_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E10_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E11_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E11_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E12_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E12_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E13_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E13_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E14_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E14_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E15_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E15_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E16_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_subset_E16_5",store.path,24,24)





### recluster the Heart, to filter the extra cell if possible.
### It doesn't work!!!!!
reduction_sample_filter <- function(sample) {
  sample  <- RunPCA(sample , assay = "SCT", verbose = F)
  sample  <- FindNeighbors(sample , reduction = "pca", dims = 1:30, verbose = F)
  sample  <- FindClusters(sample , resolution = 0.5 ,verbose = F,)
  sample  <- RunUMAP(sample , dims = 1:30, verbose = F)
  return(sample)
}


E9_5_heart <- reduction_sample_filter(E9_5_heart)
E10_5_heart <- reduction_sample_filter(E10_5_heart)
E11_5_heart <- reduction_sample_filter(E11_5_heart)
E12_5_heart <- reduction_sample_filter(E12_5_heart)
E13_5_heart <- reduction_sample_filter(E13_5_heart)
E14_5_heart <- reduction_sample_filter(E14_5_heart)
E15_5_heart <- reduction_sample_filter(E15_5_heart)
E16_5_heart <- reduction_sample_filter(E16_5_heart)


Idents(E9_5_heart) <- "SCT_snn_res.0.5"
Idents(E10_5_heart) <- "SCT_snn_res.0.5"
Idents(E11_5_heart) <- "SCT_snn_res.0.5"
Idents(E12_5_heart) <- "SCT_snn_res.0.5"
Idents(E13_5_heart) <- "SCT_snn_res.0.5"
Idents(E14_5_heart) <- "SCT_snn_res.0.5"
Idents(E15_5_heart) <- "SCT_snn_res.0.5"
Idents(E16_5_heart) <- "SCT_snn_res.0.5"

SpatialSample <- SpatialDimPlot(E9_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E9_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E10_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E10_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E11_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E11_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E12_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E12_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E13_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E13_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E14_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E14_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E15_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E15_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E16_5_heart, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E16_5",store.path,24,24)




### Keep the cell which has least 3 neighbours
### Works!!!

filter_3neighbour_cell <- function(object){
  data <- object@images$slice1@coordinates[4:5]
  dist_euclidean <- dist(data,method="euclidean")
  dist_matrix <- as.matrix(dist_euclidean)
  dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]
  counts = matrix(ncol = ncol(dist_matrix), nrow =1)
  for (i in 1:ncol(dist_matrix)){
    counts[i] <- sum(dist_matrix[,i] < 2.2)
    }
  colnames(counts) <- colnames(dist_matrix)

  filter_object <- counts[,counts>3]
  filter_barcodes <- as.character(names(filter_object))
  mainpart <- subset(object, cells = filter_barcodes)
  return(mainpart)
  }


E9_5_filter <- filter_3neighbour_cell(E9_5_heart)
E10_5_filter <- filter_3neighbour_cell(E10_5_heart)
E11_5_filter <- filter_3neighbour_cell(E11_5_heart)
E12_5_filter <- filter_3neighbour_cell(E12_5_heart)
E13_5_filter <- filter_3neighbour_cell(E13_5_heart)
E14_5_filter <- filter_3neighbour_cell(E14_5_heart)
E15_5_filter <- filter_3neighbour_cell(E15_5_heart)
E16_5_filter <- filter_3neighbour_cell(E16_5_heart)



SpatialSample <- SpatialDimPlot(E9_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E9_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E10_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E10_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E11_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E11_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E12_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E12_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E13_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E13_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E14_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E14_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E15_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E15_5",store.path,24,24)

SpatialSample <- SpatialDimPlot(E16_5_filter, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"Spatialplot_Heart_filter_E16_5",store.path,24,24)


saveRDS(E9_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5_heart_filter.rds")
saveRDS(E10_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5_heart_filter.rds")
saveRDS(E11_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5_heart_filter.rds")
saveRDS(E12_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5_heart_filter.rds")
saveRDS(E13_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E13_5_heart_filter.rds")
saveRDS(E14_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E14_5_heart_filter.rds")
saveRDS(E15_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E15_5_heart_filter.rds")
saveRDS(E16_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E16_5_heart_filter.rds")




#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### All sample integration

Mouse_embryo_all_stage.list <- list(E9_5_heart,E10_5_heart,E11_5_heart,E12_5_heart,E13_5_heart,E14_5_heart,E15_5_heart,E16_5_heart)

genes.common <- Reduce(intersect, list(rownames(E9_5_heart), rownames(E10_5_heart), rownames(E11_5_heart), rownames(E12_5_heart), rownames(E13_5_heart), rownames(E14_5_heart), rownames(E15_5_heart), rownames(E16_5_heart)))

list.features <- SelectIntegrationFeatures(object.list = Mouse_embryo_all_stage.list,
                                           nfeatures = 3000,
                                           assay = c("SCT", "SCT", "SCT", "SCT","SCT","SCT","SCT","SCT"))

Mouse_embryo_all_stage.list <- PrepSCTIntegration(object.list = Mouse_embryo_all_stage.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = F)
Mouse_embryo_all_stage.anchors <- FindIntegrationAnchors(object.list = Mouse_embryo_all_stage.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = F)                       

saveRDS(Mouse_embryo_all_stage.anchors,"/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated_anchors.rds")
saveRDS(genes.common,"/biomedja01/disk1/lmx_home/thesis/processed_data/all_stage_genes.common.rds")

Mouse_embryo_all_stage.anchors <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated_anchors.rds")
genes.common <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/all_stage_genes.common.rds")

Mouse_embryo_all_stage <- IntegrateData(anchorset = Mouse_embryo_all_stage.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = T)


# ### Xomics to process the data
# Mouse_embryo_all_stage.anchors <- readRDS("/groups/cgsd/andyliang/thesis/Heart_all_stage_integrated_anchors.rds")
# genes.common <- readRDS("/groups/cgsd/andyliang/thesis/genes.common.rds")



# saveRDS(Mouse_embryo_all_stage,"/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated.rds")
# Mouse_embryo_all_stage <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated.rds")


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

saveRDS(E9_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5_heart_filter.rds")
saveRDS(E10_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5_heart_filter.rds")
saveRDS(E11_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5_heart_filter.rds")
saveRDS(E12_5_filter,"/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5_heart_filter.rds")

E9_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E9_5_heart_filter.rds")
E10_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E10_5_heart_filter.rds")
E11_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E11_5_heart_filter.rds")
E12_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E12_5_heart_filter.rds")

### 4 stage level

Mouse_embryo_all_stage.list <- list(E9_5_heart,E10_5_heart,E11_5_heart,E12_5_heart)

genes.common <- Reduce(intersect, list(rownames(E9_5_heart), rownames(E10_5_heart), rownames(E11_5_heart), rownames(E12_5_heart)))

list.features <- SelectIntegrationFeatures(object.list = Mouse_embryo_all_stage.list,
                                           nfeatures = 3000,
                                           assay = c("SCT", "SCT", "SCT", "SCT"))

Mouse_embryo_all_stage.list <- PrepSCTIntegration(object.list = Mouse_embryo_all_stage.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = T)
Mouse_embryo_all_stage.anchors <- FindIntegrationAnchors(object.list = Mouse_embryo_all_stage.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = T)           

Mouse_embryo_all_stage <- IntegrateData(anchorset = Mouse_embryo_all_stage.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = T)


saveRDS(Mouse_embryo_all_stage,"/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated_9_12.rds")



### Switch to integrated assays
Mouse_embryo_all_stage <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated_9_12.rds")
DefaultAssay(Mouse_embryo_all_stage) <- "integrated"

# Rerun dimensionality reduction and clustering on integrated object.
Mouse_embryo_all_stage <- RunPCA(Mouse_embryo_all_stage, verbose = F)
Mouse_embryo_all_stage <- FindNeighbors(Mouse_embryo_all_stage, reduction = "pca", dims = 1:30)
Mouse_embryo_all_stage <- FindClusters(Mouse_embryo_all_stage, resolution = 0.5, verbose = F)
Mouse_embryo_all_stage <- RunUMAP(Mouse_embryo_all_stage, reduction = "pca", dims = 1:30)


DefaultAssay(Mouse_embryo_all_stage) <- "SCT"
Idents(Mouse_embryo_all_stage) <- "timepoint"
p1 <- DimPlot(Mouse_embryo_all_stage, reduction = "umap", label = TRUE)
saveplot(plot_name,p1,"umap_timepoint",store.path,5,5)


Idents(Mouse_embryo_all_stage) <- "integrated_snn_res.0.5"
p1 <- DimPlot(Mouse_embryo_all_stage, reduction = "umap", label = TRUE)
saveplot(plot_name,p1,"umap_res0.5",store.path,5,5)

Mouse_embryo_all_stage <- RunTSNE(Mouse_embryo_all_stage, reduction = "pca", dims = 1:30)
p1 <- DimPlot(Mouse_embryo_all_stage, reduction = "tsne", label =TRUE)
saveplot(plot_name,p1,"tsne_res0.5",store.path,5,5)




# Find the cluster specfic gene

Mouse_embryo_all_stage <- PrepSCTFindMarkers(Mouse_embryo_all_stage, assay = "SCT", verbose = TRUE)### before FindMarker, need to prepSCT
saveRDS(Mouse_embryo_all_stage,"/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated_9_12_prepSCT.rds")
Mouse_embryo_all_stage <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_all_stage_integrated_9_12_prepSCT.rds")

Idents(Mouse_embryo_all_stage) <- "integrated_snn_res.0.5"
Mouse_embryo_all_stage.markers <- FindAllMarkers(Mouse_embryo_all_stage,assays="SCT" ,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,recorrect_umi = FALSE)
write.csv(Mouse_embryo_all_stage.markers,"/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_9_12_deg.csv")
Mouse_embryo_all_stage.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
p <- DoHeatmap(Mouse_embryo_all_stage, features = top10$gene) + NoLegend()
saveplot(plot_name,p,"res0.5_cluster_specfic_gene",store.path,width = 5,height = 10)




Mouse_embryo_all_stage.markers_Spatial <- FindAllMarkers(Mouse_embryo_all_stage,assays="Spatial" ,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,recorrect_umi = FALSE)
write.csv(Mouse_embryo_all_stage.markers_Spatial,"/biomedja01/disk1/lmx_home/thesis/processed_data/Heart_9_12_deg_Spatial.csv")
Mouse_embryo_all_stage.markers_Spatial %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DefaultAssay(Mouse_embryo_all_stage) <- "Spatial"
p <- DoHeatmap(Mouse_embryo_all_stage, features = top10$gene) + NoLegend()
saveplot(plot_name,p,"res0.5_Spatial_cluster_specfic_gene",store.path,width = 5,height = 10)


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

### Cell type annotation
### Based on the paper,
### Feng, W., Bais, A., He, H. et al. 
### Single-cell transcriptomic analysis identifies murine heart molecular features at embryonic and neonatal stages. 
### Nat Commun 13, 7960 (2022). https://doi.org/10.1038/s41467-022-35691-7


Dotplot_marker <- function(object,markerlist,name,width,height){
  

  p_all_markers = DotPlot(object, 
                      features = markerlist,
                      scale = T,assay='SCT',group.by = "integrated_snn_res.0.5" )+
    theme_bw()+
    scale_color_continuous(low="grey",high =  "red")+
    theme(legend.position = "right",legend.box = "vertical",
    legend.margin=margin(t= 0, unit='cm'),
    #legend.spacing = unit(0,"in"),
    axis.text.x  = element_text(color="black",size=12,angle = 45,vjust = 0.5, hjust=0.5),
    axis.text.y  = element_text(color="black",size=12),
    legend.text = element_text(size =12,color="black"),
    legend.title = element_text(size =12,color="black")
    ) 

    saveplot(plot_name,p_all_markers,name,store.path,width,height)
}


pan_CM = c("Ttn")
Atrial_CM = c("Sln")
Ventricular_CM = c("Myl2")
pan_EC = c("Pecam1")
Endo_EC = c("Npr3")
Vas_EC = c("Fab4")
Epi = c("Wt1","Tbx18","Aldh1a2")
Fb = c("Postn")
Mural = c("Pdgfrb")
Immune = c("C1q2")
Blood = c("Hb1-a1")


markerlist =list(
  pan_CM = pan_CM,
  Atrial_CM = Atrial_CM,
  Ventricular_CM = Ventricular_CM,
  pan_EC = pan_EC,
  Endo_EC = Endo_EC,
  Vas_EC = Vas_EC,
  Epi = Epi,
  Fb = Fb,
  Mural = Mural,
  Immune = Immune,
  Blood = Blood
)

Dotplot_marker(Mouse_embryo_all_stage,markerlist,"Heart_Dotplot_Mouse_embryo_all_stage_9_12",12,8)




Fb = c("Postn","Ptn","Col3a1","Nnat")
Atrial_CM = c("Myl7","Myl1","Sln","Myl4")
Ventricular_CM = c("Myl2","Myl3","Pln","Actc1","Tnnt2")
Progenitor = c("Igfbp2")
Endo = c("Aap1","Gja4","Ramp2","Egfl7","Fn1","Fabp5")
# Epi = c("Ccl25","Krt8","Krt18","Crapb2")
Erythroid = c("Hba-x","Hbb-y","Hba-a1","Hba-a2")
# lymphocyte = c("Ccl4","Napsa","Plac8","Cd74")
# Mast = c("Cpa3","Cma1")
Marcophahge = c("Gpx1","Mt1")
Mesotheial = c("Upk3b","Itm2a","Mgst1")
Myofib = c("Eln","Fbln5","Acta2","Cxcl12")

markerlist =list(
  Fb = Fb,
  Atrial_CM = Atrial_CM,
  Ventricular_CM = Ventricular_CM,
  Progenitor = Progenitor,
  Endo = Endo,
  # Epi = Epi,
  Erythroid = Erythroid,
  # lymphocyte = lymphocyte,
  # Mast = Mast,
  Marcophahge = Marcophahge,
  Mesotheial = Mesotheial,
  Myofib = Myofib
)

Dotplot_marker(Mouse_embryo_all_stage,markerlist,"Heart_Dotplot_Mouse_embryo_all_stage_9_12_MCA",18,8)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



### Mouse_embryo_all_stage annotation

Mouse_embryo_all_stage@meta.data$Cell_Type <- Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5


Mouse_embryo_all_stage@meta.data$Cell_Type <- "Other"


Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="0",]$Cell_Type <- "Ventricular_CM"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="1",]$Cell_Type <- "Fibroblast"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="2",]$Cell_Type <- "Atrial_CM"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="3",]$Cell_Type <- "Ventricular_CM"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="4",]$Cell_Type <- "Blood"

# Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="5",]$Cell_Type <- "Fb_Postn_high"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="5",]$Cell_Type <- "Fibroblast"

# Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="6",]$Cell_Type <- "Fb_Ptn_high"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="6",]$Cell_Type <- "Fibroblast"

Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="7",]$Cell_Type <- "Fibroblast"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="8",]$Cell_Type <- "Ventricular_CM"
# Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="9",]$Cell_Type <- "Fb_Alnat_high"
Mouse_embryo_all_stage@meta.data[Mouse_embryo_all_stage@meta.data$integrated_snn_res.0.5=="9",]$Cell_Type <- "Fibroblast"



### draw the new umap with circle
tSNE<- as.data.frame(Mouse_embryo_all_stage@reductions$umap@cell.embeddings)
Idents(Mouse_embryo_all_stage) <-"Cell_Type"
Cluster <- Idents(Mouse_embryo_all_stage)
table(Cluster)
tSNE <- cbind(tSNE,Cluster)
#建立自定义主题：
mytheme <- theme_void() + #空白主题，便于我们后期添加tSNE箭头
  theme(plot.margin = margin(5.5,15,5.5,5.5)) #画布空白页缘调整

#建立映射，添加散点：
p <- ggplot(data = tSNE, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster),
             size = 0.4,
             alpha = 0.8)
p

p1 <- p +
  stat_ellipse(aes(color = Cluster),
               level = 0.95, linetype = 1, show.legend = F) +
  mytheme

p2 <- p +
  stat_ellipse(aes(color = Cluster),
               level = 0.95, linetype = 2, show.legend = F) +
  mytheme


p3 <- p +
  stat_ellipse(aes(color = Cluster, fill = Cluster),
               level = 0.95, linetype = 1, show.legend = F,
               geom = 'polygon', alpha = 0.1) +
  mytheme

p4 <- p3 +
  theme_dr(xlength = 0.2, #x轴长度
           ylength = 0.2, #y轴长度
           arrow = grid::arrow(length = unit(0.1, "inches"), #箭头大小/长度
                               ends = 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())

label <- tSNE %>%
  group_by(Cluster)%>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

p5 <- p4 +
  geom_text(data = label,
            aes(x = UMAP_1, y = UMAP_2, label = Cluster,
            fontface = "bold", #粗体强调
            color = 'black', size = 4))

p6 <- p5 +
  guides(color = guide_legend(override.aes = list(size = 5)))


saveplot(plot_name,p6,"umap_new",store.path,8,5)


### Observe the cell type localization in spatial
Idents(Mouse_embryo_all_stage) <-"Cell_Type"
Mouse_embryo_all_stage.list <- SplitObject(Mouse_embryo_all_stage, split.by = "timepoint")
Idents(Mouse_embryo_all_stage) <-"Cell_Type"
SpatialSample <- SpatialDimPlot(Mouse_embryo_all_stage, crop = FALSE, label = FALSE, pt.size.factor = 1.5, label.size = 4, ncol = 2)
saveplot(plot_name,SpatialSample,"Spatialplot_cellannotation",store.path,10,10)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### After cell type annotation of it, observe the cell type composition of each sample.
### Bar plot to show

# data<- data.frame(Cell_Type=rep(c("pan_CM","Atrial_CM","Ventricular_CM","pan_EC","Endo_EC","Vas_EC","Epi","Fb","Mural","Immune","Blood"),8),
#   sample_group = rep(c("Control","VclcKO"),each = 4),
#   count = c(0.077633592,0.5010258,0.421052632,0.000287976,0.203794773,0.422249852,0.35483871,0.019116665))

Mouse_embryo_all_stage.list <- SplitObject(Mouse_embryo_all_stage, split.by = "timepoint")

Composistion_E9_5 <- table(Mouse_embryo_all_stage.list$E9.5@meta.data$Cell_Type)/sum(table(Mouse_embryo_all_stage.list$E9.5@meta.data$Cell_Type))
Composistion_E10_5 <- table(Mouse_embryo_all_stage.list$E10.5@meta.data$Cell_Type)/sum(table(Mouse_embryo_all_stage.list$E10.5@meta.data$Cell_Type))
Composistion_E11_5 <- table(Mouse_embryo_all_stage.list$E11.5@meta.data$Cell_Type)/sum(table(Mouse_embryo_all_stage.list$E11.5@meta.data$Cell_Type))
Composistion_E12_5 <- table(Mouse_embryo_all_stage.list$E12.5@meta.data$Cell_Type)/sum(table(Mouse_embryo_all_stage.list$E12.5@meta.data$Cell_Type))




data<- data.frame(Cell_Type=rep(c("Ventricular_CM","Atrial_CM","Fibroblast","Blood"),4),
  sample_group = rep(c("E9_5","E10_5","E11_5","E12_5"),each = 4),
  count = c(Composistion_E9_5["Ventricular_CM"],Composistion_E9_5["Atrial_CM"],Composistion_E9_5["Fibroblast"],Composistion_E9_5["Blood"],
Composistion_E10_5["Ventricular_CM"],Composistion_E10_5["Atrial_CM"],Composistion_E10_5["Fibroblast"],Composistion_E10_5["Blood"],
Composistion_E11_5["Ventricular_CM"],Composistion_E11_5["Atrial_CM"],Composistion_E11_5["Fibroblast"],Composistion_E11_5["Blood"],
Composistion_E12_5["Ventricular_CM"],Composistion_E12_5["Atrial_CM"],Composistion_E12_5["Fibroblast"],Composistion_E12_5["Blood"]
))


data$Cell_Type <- factor(data$Cell_Type, levels = c("Ventricular_CM","Atrial_CM","Fibroblast","Blood"))


gg <- ggplot(data, aes(x = sample_group, y = count, fill = Cell_Type))+
  geom_bar(stat = "identity", width = 0.7)+ # 柱状图绘制
  geom_flow(aes(alluvium = Cell_Type), alpha = 0.5) + # 添加柱状图后的条带
  labs(fill = "Cell Type")+
  theme_bw()+ # 将主题调整为白色背景和浅灰色网格线
  xlab("")+ # 去掉x轴的标题
  ylab("Relative Cell Composition (%)")+ # 设置y轴的标签
  scale_x_discrete(labels = c("E9_5","E10_5","E11_5","E12_5"))+
  theme(panel.grid.major.x = element_blank(),
      axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size=rel(0.8)),
        legend.text = element_text(size = rel(1),hjust =0))+ # 更改x轴、y轴的字体大小、刻度线等
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0.05))+
  geom_text(aes(label = paste0(round(count*100), "%")), position = position_stack(vjust = 0.5), size = 3)


saveplot(plot_name,gg,"percentage",store.path,10,10)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### check speficifc gene expression



DefaultAssay(Mouse_embryo_all_stage) <- "SCT"
Sample <- SpatialFeaturePlot(object = Mouse_embryo_all_stage, features = c("Nkx2-5","Bnc1","Nr2f2","Hey2","Hand2","Mesp1","Hopx","Prdm16"), ncol = 4,pt.size.factor = 4) #alpha = c(0.1, 1)
saveplot(plot_name,Sample,"SpatialFeaturePlot_",store.path,24,48)



###Mesp1 was critical in regulating the expression of several key genes during heart development, including the EMT genes Snai1 and Zeb2, 
### the migration gene Rasgrp3, and the cardiac cell commitment genes Etv2, Hand1, Myl7, Gata4, Flk1, and Pdgfra
DefaultAssay(Mouse_embryo_all_stage) <- "SCT"
Sample <- SpatialFeaturePlot(object = Mouse_embryo_all_stage, features = c("Snail1","Zeb2","Rasgrp3","Etv2","Hand1","Myl7","Gata4","Flk1","Pdgfra"), ncol = 4,pt.size.factor = 4) #alpha = c(0.1, 1)
saveplot(plot_name,Sample,"SpatialFeaturePlot_Mesp1_related",store.path,24,48)


###Prdm16 cardiomyocyte-specific knockout mouse heart and identified Pdrm16 as an important transcription factor in myocardial densification.
###the transcriptional regulation of Prdm16 has been found to be chamber-specific and to be associated with the enrichment of Tbx5 and Hand1 in the left ventricular region
DefaultAssay(Mouse_embryo_all_stage) <- "SCT"
Sample <- SpatialFeaturePlot(object = Mouse_embryo_all_stage, features = c("Tbx5","Hand1"), ncol = 4,pt.size.factor = 4) #alpha = c(0.1, 1)
saveplot(plot_name,Sample,"SpatialFeaturePlot_Pdrm16_related",store.path,24,24)


###Ap-2
DefaultAssay(Mouse_embryo_all_stage) <- "SCT"
Sample <- SpatialFeaturePlot(object = Mouse_embryo_all_stage, features = c("Actg2"), ncol = 4,pt.size.factor = 4) #alpha = c(0.1, 1)
saveplot(plot_name,Sample,"SpatialFeaturePlot_Cellmarker_related",store.path,24,24)





### Nkx2-5 
### Nkx2–5 was necessary for cardiomyocyte maturation and the establishment of ventricular structure
### Nkx2–5 was found to directly bind the Cxcr2 and Cxcr4 genomic loci and activate their transcription in the second heart field, 
### thereby participating in the migration and differentiation of cardiac progenitor cells

### Hand2
### Hand2 has been found to be a specific transcription factor for OFT cells (Fig. 3). In a Hand2-null mouse embryonic model, 
### OFT cardiomyocytes lost their specification properties, 
### and the differentiation and migration of right ventricular cardiomyocytes were also affected



color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))


genes_to_check = c("Nkx2-5","Bnc1","Nr2f2","Hey2","Hand2","Mesp1","Hopx","Prdm16" )
genes_to_check2<-as.data.frame(genes_to_check)

#genes_to_check2<-t(genes_to_check2)

Mouse_embryo_all_stage$timepoint <- factor(x = Mouse_embryo_all_stage$timepoint, levels = c('E9.5','E10.5','E11.5','E12.5'))

p1 <- VlnPlot(Mouse_embryo_all_stage,features =  genes_to_check,
              group.by  = "Cell_Type",
              split.by = "timepoint",
              flip = T,stack = T,cols = color
              )
p2<-p1 + NoLegend()


saveplot(plot_name,p1,"VlnPlot_genes",store.path,12,12)





#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Re-cluster for specfic cell type





### scRNA-seq validation






#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### In this part, we will zoom in the co-localization of the cell type in each stage









#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### In this part, zoom in the latest stage of heart E16.5

E16_5_heart <- readRDS("/biomedja01/disk1/lmx_home/thesis/processed_data/SCT_E16_5_heart_filter.rds")




Dotplot_marker <- function(object,markerlist,name,width,height){
  

  p_all_markers = DotPlot(object, 
                      features = markerlist,
                      scale = T,assay='SCT',group.by = "SCT_snn_res.0.5" )+
    theme_bw()+
    scale_color_continuous(low="grey",high =  "red")+
    theme(legend.position = "right",legend.box = "vertical",
    legend.margin=margin(t= 0, unit='cm'),
    #legend.spacing = unit(0,"in"),
    axis.text.x  = element_text(color="black",size=12,angle = 45,vjust = 0.5, hjust=0.5),
    axis.text.y  = element_text(color="black",size=12),
    legend.text = element_text(size =12,color="black"),
    legend.title = element_text(size =12,color="black")
    ) 

    saveplot(plot_name,p_all_markers,name,store.path,width,height)
}



Fb = c("Postn","Ptn","Col3a1","Nnat")
Atrial_CM = c("Myl7","Myl1","Sln","Myl4")
Ventricular_CM = c("Myl2","Myl3","Pln","Actc1","Tnnt2")
Progenitor = c("Igfbp2")
Endo = c("Gja4","Ramp2","Egfl7","Fn1","Fabp5")
Epi = c("Ccl25","Krt8","Krt18","Crapb2")
Erythroid = c("Hba-x","Hbb-y","Hba-a1","Hba-a2")
# lymphocyte = c("Ccl4","Napsa","Plac8","Cd74")
# Mast = c("Cpa3","Cma1")
Marcophahge = c("Gpx1","Mt1")
Mesotheial = c("Upk3b","Itm2a","Mgst1")
Myofib = c("Eln","Fbln5","Acta2","Cxcl12")

markerlist =list(
  Fb = Fb,
  Atrial_CM = Atrial_CM,
  Ventricular_CM = Ventricular_CM,
  Progenitor = Progenitor,
  Endo = Endo,
  Epi = Epi,
  Erythroid = Erythroid,
  # lymphocyte = lymphocyte,
  # Mast = Mast,
  Marcophahge = Marcophahge,
  Mesotheial = Mesotheial,
  Myofib = Myofib
)

Dotplot_marker(E16_5_heart,markerlist,"Heart_Dotplot_E16_5_MCA",20,8)































