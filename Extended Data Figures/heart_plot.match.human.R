## Plot the UMAP of the dataset
library(Seurat)
library(dplyr)
inputdir <- "./load_files/"
outputdir <- "./report/"



##https://www.nature.com/articles/s41586-020-2797-4#data-availability
human <- readRDS(file = "./human_heart/global_raw_seurat.rds")
exp_human <- rownames(human)[Matrix::rowSums(human$RNA@data!=0) >= 20]
## Normalize data [The data has not been normalized]
human <- NormalizeData(human)


Idents(human) <- "cell_type"
avgs_human <- log(AverageExpression(human)$RNA + 1)
cls_order <- c("Atrial_Cardiomyocyte","Ventricular_Cardiomyocyte","Myeloid","Lymphoid","Mesothelial","Endothelial","Pericytes","Smooth_muscle_cells","Fibroblast","Adipocytes","Neuronal")



pig <- readRDS(file = paste0(inputdir, "OE.heart.all.seurat.filter.final.rds"))
hvg_pig <- pig[, as.character(pig@meta.data$condition) == 'h0'] %>%
            SplitObject(., split.by = "samplename") %>%
            lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
            SelectIntegrationFeatures(., nfeatures = 2000)

Idents(pig) <- "cluster"
avg_pig <- log(AverageExpression(pig)$RNA + 1) 
##exp_pig <- rownames(pig)[Matrix::rowSums(pig$RNA@counts!=0) >= 15]


pig_order <- c("Cardiomyocytes","Myeloid","Lymphoid","Mesothelial","Endothelial","Pericytes","Smooth Muscle","Fibroblasts","Neuronal")


gene_use <- intersect(hvg_pig, exp_human)
cor_mat <- cor(avg_pig[gene_use, ], avgs_human[gene_use, ], method = "p")


pdf(paste0(outputdir, "SF1.Heart.human-pig.heatmap.pdf"), width = 7, height = 5)
pheatmap::pheatmap(cor_mat[pig_order, cls_order], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis::viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 5, fontsize_row = 8)
dev.off()



##------------------------------------------------------------------------
## UMAP visualization
source("~/project/HIP/scripts/hip.fun.R")
pig <- readRDS(file = paste0(inputdir, "OE.heart.all.seurat.filter.final.rds"))

pig@meta.data$group <- gsub("VLMC", "Fibroblasts", pig@meta.data$group)
gp_colors <- c("#4472C4", "#70AD47", "#E0BC8D", "#FFC000", "#C83200") %>% 
          setNames(., c("Cardiomyocytes", "Neuronal", "Fibroblasts", "Endo_mural", "immune"))
DimFig(pig, group.by = "group", cols = gp_colors, file_name = "SF1.Heart.lowreso.UMAP", plot.scale = 0.9)




pig@meta.data$cluster <- gsub("Endothelial", "Endothelial-2", pig@meta.data$cluster) %>%
				gsub("Mesothelial", "Endothelial-1", .)
cls_cols <- setNames(c("#6650A2", "#892501", "#d15957", "#ea7807", "#f9c159", "#db9130", "#ce871c", "#e8bafc", "#e1bd8e"), c("Cardiomyocytes","Myeloid","Lymphoid","Endothelial-1","Endothelial-2","Pericytes","Smooth Muscle","Fibroblasts","Neuronal"))
DimFig(pig, group.by = "cluster", cols = cls_cols, file_name = "SF1.Heart.highreso.UMAP", plot.scale = 0.9)














