## Plot the UMAP of the dataset
source("~/project/HIP/scripts/hip.fun.R")
library(ggdendro)


##https://www.science.org/doi/10.1126/science.aat5031
human <- readRDS(file = "./human_kidney/human_kidney_0504_2021.rds")
exp_human <- rownames(human)[Matrix::rowSums(human$RNA@data!=0) >= 20]


Idents(human) <- "broad_celltype"
avgs_human <- log(AverageExpression(human)$RNA + 1)
cls_order <- c("Proximal tubule","Connecting tubule","Principal cell","Thick ascending limb of Loop of Henle","Intercalated cell","Transitional urothelium","Pelvic epithelium","Epithelial progenitor cell","Podocyte","MNP-a/classical monocyte derived", "MNP-b/non-classical monocyte derived", "MNP-c/dendritic cell", "MNP-d/Tissue macrophage","Neutrophil","Plasmacytoid dendritic cell","Mast cell","B cell","CD4 T cell","CD8 T cell","NK cell","NKT cell","Peritubular capillary endothelium","Ascending vasa recta endothelium","Glomerular endothelium","Descending vasa recta endothelium","Fibroblast","Myofibroblast")    

pig <- readRDS(file = paste0(inputdir, "OE.kidney.all.seurat.filter.final.rds"))
hvg_pig <- subset(pig, condition == 'h0') %>%
            SplitObject(., split.by = "samplename") %>%
            lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
            SelectIntegrationFeatures(., nfeatures = 2000)

Idents(pig) <- "cluster"
avg_pig <- log(AverageExpression(pig)$RNA + 1)
exp_pig <- rownames(pig)[Matrix::rowSums(pig$RNA@counts!=0) >= 15]


pig_order <- c("Proximal Tubule","CNT PC","Loop of Henle","Intercalated Cell","Pedocyte","Myeloid","Lymphoid","Endothelium","Fibroblast")


gene_use <- intersect(hvg_pig, exp_human)
cor_mat <- cor(avg_pig[gene_use, ], avgs_human[gene_use, ], method = "p")


pdf(paste0(outputdir, "SF1.Kidney.human-pig.heatmap.pdf"), width = 7, height = 5)
pheatmap::pheatmap(cor_mat[pig_order, cls_order], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 5, fontsize_row = 8)
dev.off()






##------------------------------------------------------------------------
## UMAP visualization
source("~/project/HIP/scripts/hip.fun.R")
pig <- readRDS(file = paste0(inputdir, "OE.kidney.all.seurat.filter.final.rds"))
gp_colors <- c("#4472C4", "#70AD47", "#E0BC8D", "#FFC000")%>% 
          setNames(., c(c("Nephron", "Immune", "Endothelium", "Stroma")))
DimFig(pig, group.by = "group", cols = gp_colors, file_name = "SF1.Kidney.lowreso.UMAP", plot.scale = 0.9)



cls_cols <- setNames(c("#6016ff", "#5d00aa", "#842bc4", "#4a08a0", "#5118b2", "#fff187", "#ffffa8", "#DFD4E9", "#FFA500"), 
	c("Proximal Tubule","CNT PC","Loop of Henle","Intercalated Cell","Pedocyte","Myeloid","Lymphoid","Endothelium","Fibroblast"))
DimFig(pig, group.by = "cluster", cols = cls_cols, file_name = "SF1.Kidney.highreso.UMAP", plot.scale = 0.9)





