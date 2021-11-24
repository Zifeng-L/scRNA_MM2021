# Primary analysis of all cells from 18 samples using seurat.

rm(list =ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
options(future.globals.maxSize = 1000000 * 1024^2)
data_1 <- Read10X(data.dir = "../MM01_D/gz")
data_2 <- Read10X(data.dir = "../MM02_D/gz")
data_3 <- Read10X(data.dir = "../MM03_D/gz")
data_4 <- Read10X(data.dir = "../MM04_D/gz")
data_5 <- Read10X(data.dir = "../MM05_D/gz")
data_6 <- Read10X(data.dir = "../MM06_D/gz")
data_7 <- Read10X(data.dir = "../MM07_D/gz")
data_8 <- Read10X(data.dir = "../MM08_D/gz")
data_9 <- Read10X(data.dir = "../MM09_D/gz")
data_10 <- Read10X(data.dir = "../MM10_D/gz")
data_11 <- Read10X(data.dir = "../MM11_D/gz")
data_12 <- Read10X(data.dir = "../MM12_D/gz")

data_13 <- Read10X(data.dir = "../MM13_R/gz")
data_14 <- Read10X(data.dir = "../MM14_R/gz")
data_15 <- Read10X(data.dir = "../MM15_R/gz")
data_16 <- Read10X(data.dir = "../MM16_R/gz")
data_17 <- Read10X(data.dir = "../MM17_R/gz")
data_18 <- Read10X(data.dir = "../MM18_R/gz")

data_1 <- CreateSeuratObject(data_1, project = "MM01_D")
data_2 <- CreateSeuratObject(data_2, project = "MM02_D")
data_3 <- CreateSeuratObject(data_3, project = "MM03_D")
data_4 <- CreateSeuratObject(data_4, project = "MM04_D")
data_5 <- CreateSeuratObject(data_5, project = "MM05_D")
data_6 <- CreateSeuratObject(data_6, project = "MM06_D")
data_7 <- CreateSeuratObject(data_7, project = "MM07_D")
data_8 <- CreateSeuratObject(data_8, project = "MM08_D")
data_9 <- CreateSeuratObject(data_9, project = "MM09_D")
data_10 <- CreateSeuratObject(data_10, project = "MM10_D")
data_11 <- CreateSeuratObject(data_11, project = "MM11_D")
data_12 <- CreateSeuratObject(data_12, project = "MM12_D")

data_13 <- CreateSeuratObject(data_13, project = "MM13_R")
data_14 <- CreateSeuratObject(data_14, project = "MM14_R")
data_15 <- CreateSeuratObject(data_15, project = "MM15_R")
data_16 <- CreateSeuratObject(data_16, project = "MM16_R")
data_17 <- CreateSeuratObject(data_17, project = "MM17_R")
data_18 <- CreateSeuratObject(data_18, project = "MM18_R")

data_1@meta.data$group <- "MM01_D"
data_2@meta.data$group <- "MM02_D"
data_3@meta.data$group <- "MM03_D"
data_4@meta.data$group <- "MM04_D"
data_5@meta.data$group <- "MM05_D"
data_6@meta.data$group <- "MM06_D"
data_7@meta.data$group <- "MM07_D"
data_8@meta.data$group <- "MM08_D"
data_9@meta.data$group <- "MM09_D"
data_10@meta.data$group <- "MM10_D"
data_11@meta.data$group <- "MM11_D"
data_12@meta.data$group <- "MM12_D"

data_13@meta.data$group <- "MM13_R"
data_14@meta.data$group <- "MM14_R"
data_15@meta.data$group <- "MM15_R"
data_16@meta.data$group <- "MM16_R"
data_17@meta.data$group <- "MM17_R"
data_18@meta.data$group <- "MM18_R"

data_1@meta.data$category <- "Diagnosis"
data_2@meta.data$category <- "Diagnosis"
data_3@meta.data$category <- "Diagnosis"
data_4@meta.data$category <- "Diagnosis"
data_5@meta.data$category <- "Diagnosis"
data_6@meta.data$category <- "Diagnosis"
data_7@meta.data$category <- "Diagnosis"
data_8@meta.data$category <- "Diagnosis"
data_9@meta.data$category <- "Diagnosis"
data_10@meta.data$category <- "Diagnosis"
data_11@meta.data$category <- "Diagnosis"
data_12@meta.data$category <- "Diagnosis"

data_13@meta.data$category <- "Relapse"
data_14@meta.data$category <- "Relapse"
data_15@meta.data$category <- "Relapse"
data_16@meta.data$category <- "Relapse"
data_17@meta.data$category <- "Relapse"
data_18@meta.data$category <- "Relapse"

data_merge <- merge(data_1,c(data_2,data_3,data_4,data_5,data_6,data_7,data_8, data_9, data_10,data_11, data_12,data_13,data_14,data_15,data_16,data_17,data_18),project = "merge",merge.data = TRUE)

data_list = SplitObject(data_merge, split.by = "group")

for (i in 1:length(data_list)) {
  data_list[[i]][["percent.mt"]]<-PercentageFeatureSet(data_list[[i]],pattern="^MT-")
  data_list[[i]][["HBB"]]<-PercentageFeatureSet(data_list[[i]],pattern="^HBB")
  data_list[[i]] <- subset(data_list[[i]], subset =  nFeature_RNA >500 & nFeature_RNA < 6000 & percent.mt <10 & HBB < 2 & nCount_RNA > 2000 )
  data_list[[i]] <- data_list[[i]][ ! grepl("MALAT1", rownames(data_list[[i]])), ]
  data_list[[i]] <- data_list[[i]][ ! grepl("^MT-", rownames(data_list[[i]])), ]
  data_list[[i]] <- data_list[[i]][ ! grepl("^RP[SL]", rownames(data_list[[i]])), ]
  data_list[[i]] <- data_list[[i]][ ! grepl("^HB[AB]", rownames(data_list[[i]])), ]
  data_list[[i]] <- data_list[[i]][ ! grepl("^IGKV", rownames(data_list[[i]])), ]
  data_list[[i]] <- data_list[[i]][ ! grepl("^IGLV", rownames(data_list[[i]])), ]
  data_list[[i]] <- data_list[[i]][ ! grepl("^IGHV", rownames(data_list[[i]])), ]
  data_list[[i]] <- NormalizeData(data_list[[i]], verbose = FALSE)
  data_list[[i]] <- FindVariableFeatures(data_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
data.anchors <- FindIntegrationAnchors(object.list = data_list, anchor.features = 2000,dims = 1:30)
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:30)
DefaultAssay(data.integrated) <- "integrated"
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30,resolution=0.3)
data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:30) %>% FindClusters()

saveRDS(data.integrated,file = "../immu_CCA18.rds")

