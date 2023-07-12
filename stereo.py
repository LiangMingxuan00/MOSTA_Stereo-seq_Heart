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

###saveplot(): save the pdf file
saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
	ggsave(filename, plot, path = store.path, width = width, height = height)
}
store.path <- c("/biomedja01/disk1/lmx_home/thesis/processed_data")
plot_name <- c("")

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
Dotplot_marker(E15_5,markerlist,"Heart_Dotplot_E15_5",8,4)
Dotplot_marker(E16_5,markerlist,"Heart_Dotplot_E16_5",8,4)




### observe and subset Heart cluster based on Dotplot

obs_marker <- function(object){
	object@meta.data$new_ident <- object@meta.data$SCT_snn_res.0.5
	object@meta.data$new_ident <- "Other"
	object@meta.data[object@meta.data$SCT_snn_res.0.5=="5",]$new_ident <- "Heart"
}

E9_5@meta.data$new_ident <- E9_5@meta.data$SCT_snn_res.0.5
E9_5@meta.data$new_ident <- "Other"
E9_5@meta.data[E9_5@meta.data$SCT_snn_res.0.5=="5",]$new_ident <- "Heart"

Idents(E9_5) <- "new_ident"
SpatialSample <- SpatialDimPlot(E9_5, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 1)
saveplot(plot_name,SpatialSample,"E9_5_Heart",store.path,24,24)



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
E11_5_heart <- subset(E11_5,idents=c("5"))
E12_5_heart <- subset(E12_5,idents=c("5"))
E13_5_heart <- subset(E13_5,idents=c("5"))
E14_5_heart <- subset(E14_5,idents=c("5"))
E15_5_heart <- subset(E15_5,idents=c("5"))
E16_5_heart <- subset(E16_5,idents=c("5"))



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
###R script
### In this part, we will integration all stage heart sample together to identify the cell types

Mouse_embryo_all_stage.list <- list(E9_5= E9_5_heart,E10_5=E10_5_heart,E10_5=E10_5_heart,E11_5=E11_5_heart,E12_5=E12_5_heart,E13_5=E13_5_heart,E14_5=E14_5_heart,E15_5=E15_5_heart,E16_5=E16_5_heart)

genes.common <- Reduce(intersect, list(rownames(E9_5_heart), rownames(E10_5_heart), rownames(E11_5_heart), rownames(E12_5_heart), rownames(E13_5_heart), rownames(E14_5_heart), rownames(E15_5_heart), rownames(E16_5_heart)))

list.features <- SelectIntegrationFeatures(object.list = Mouse_embryo_all_stage.list,
                                           nfeatures = 2000,
                                           assay = c("SCT", "SCT", "SCT", "SCT","SCT","SCT","SCT","SCT"))
Mouse_embryo_all_stage.list <- PrepSCTIntegration(object.list = Mouse_embryo_all_stage.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = F)
Mouse_embryo_all_stage.anchors <- FindIntegrationAnchors(object.list = Mouse_embryo_all_stage.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = F)                       
Mouse_embryo_all_stage <- IntegrateData(anchorset = Mouse_embryo_all_stage.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = F)


### Switch to integrated assays
DefaultAssay(Mouse_embryo_all_stage) <- "integrated"

# Rerun dimensionality reduction and clustering on integrated object.
Mouse_embryo_all_stage <- RunPCA(Mouse_embryo_all_stage, verbose = F)
Mouse_embryo_all_stage <- FindNeighbors(Mouse_embryo_all_stage, reduction = "pca", dims = 1:30)
Mouse_embryo_all_stage <- FindClusters(Mouse_embryo_all_stage, resolution = 0.5, verbose = F)
Mouse_embryo_all_stage <- RunUMAP(Mouse_embryo_all_stage, reduction = "pca", dims = 1:30)





### Cell type annotation
markerlist <- c("")
Dotplot_marker(Mouse_embryo_all_stage,markerlist,"Heart_Dotplot_Mouse_embryo_all_stage")



### After cell type annotation of it, observe the cell type composition of each sample.
### Bar plot to show

data<- data.frame(Cell_Type=rep(c("MESCNCC","VSMCCNCC","SHF","Other"),2),
	sample_group = rep(c("Control","VclcKO"),each = 4),
	count = c(0.077633592,0.5010258,0.421052632,0.000287976,0.203794773,0.422249852,0.35483871,0.019116665))

data$Cell_Type <- factor(data$Cell_Type, levels = c("VSMCCNCC", "MESCNCC", "SHF", "Other"))


gg <- ggplot(data, aes(x = sample_group, y = count, fill = Cell_Type))+
  geom_bar(stat = "identity", width = 0.7)+ # 柱状图绘制
  geom_flow(aes(alluvium = Cell_Type), alpha = 0.5) + # 添加柱状图后的条带
  scale_fill_manual(values = c("SHF"="#FE3CDB","VSMCCNCC"="#548235","MESCNCC"="#C5E0B4","Other"="#D9D9D9"),labels = c(expression(VSMC[CNCC]),expression(MES[CNCC]),"SHF","Other"))+
  labs(fill = "Cell Type")+
  theme_bw()+ # 将主题调整为白色背景和浅灰色网格线
  xlab("")+ # 去掉x轴的标题
  ylab("Relative Cell Composition (%)")+ # 设置y轴的标签
  scale_x_discrete(labels = c("Control",italic("Vcl")~"cKO"))+
  theme(panel.grid.major.x = element_blank(),
  		axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size=rel(0.8)),
        legend.text = element_text(size = rel(1),hjust =0))+ # 更改x轴、y轴的字体大小、刻度线等
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0.05))+
  geom_text(aes(label = paste0(round(count*100), "%")), position = position_stack(vjust = 0.5), size = 3)

saveplot(plot_name,gg,"percentage",store.path,4,4)





### Re-cluster for specfic cell type





### scRNA-seq validation








#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### In this part, we will zoom in the co-localization of the cell type in each stage









#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
















































