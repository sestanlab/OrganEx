source("~/project/HIP/scripts/hip.fun.R")
library(ggdendro)


##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469
mouse <- readRDS(file = "./human_liver/mouse_liver_1201_2021.rds")
exp_mouse <- rownames(mouse)[Matrix::rowSums(mouse$RNA@data!=0) >= 15]


## Get average expression [the default values are not optimal for AverageExpression function]
GetAvg <- function(object, group.by) {
    ctd <- object$RNA@data
    ident <- setNames(as.character(object@meta.data[, group.by]), colnames(object))
    all_cls <- levels(as.factor(ident))
    exps <- lapply(all_cls, function(x) {
        sub_ctd <- ctd[, ident == x, drop = FALSE]
        avg <- Matrix::rowMeans(sub_ctd)
        avg
        }) %>%
            setNames(., all_cls) %>%
            as.data.frame(., check.names = FALSE) %>%
            as.matrix()
    return(exps)
}



Idents(mouse) <- "CellType"
avgs_mouse <- log(AverageExpression(mouse)$RNA + 1)
cls_order <- c("Hepatocyte_1","Hepatocyte_2", "Hepatocyte_3", "Hepatocyte_4", "Hepatocyte_5", "Hepatocyte_6","Hepatic_Stellate_Cells","Cholangiocytes","Erythroid_Cells","Inflammatory_Macrophage","Non-inflammatory_Macrophage", "Plasma_Cells", "Mature_B_Cells","alpha-beta_T_Cells","gamma-delta_T_Cells_1", "gamma-delta_T_Cells_2", "NK-like_Cells","Central_venous_LSECs", "Periportal_LSECs", "Portal_endothelial_Cells")

pig <- readRDS(file = paste0(inputdir, "OE.liver.all.seurat.filter.final.rds"))
hvg_pig <- subset(pig, condition == 'h0') %>%
            SplitObject(., split.by = "samplename") %>%
            lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
            SelectIntegrationFeatures(., nfeatures = 2000)


pig@meta.data$newcls <- get_ident(input_ident = pig, ident_col = "seurat_clusters", file_path = paste0(inputdir, "pig.liver.hres.anno.txt"), sample_name = "pig", label_names = NULL)
Idents(pig) <- "newcls"
avg_pig <- log(AverageExpression(pig)$RNA + 1)


pig_order <- c("Hepatocytes", "Stellate", "Cholangiocytes", "myeloid", "B cells", "Plasma cells", "NK_NKT_T", "Endothelial")

gene_use <- intersect(hvg_pig, exp_mouse)
cor_mat <- cor(avg_pig[gene_use, ], avgs_mouse[gene_use, ], method = "p")


pdf(paste0(outputdir, "SF1.liver.mouse-pig.heatmap.pdf"), width = 7, height = 5)
pheatmap::pheatmap(cor_mat[pig_order, cls_order], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 5, fontsize_row = 8)
dev.off()




##------------------------------------------------------------------------------------
## UMAP visualization
source("~/project/HIP/scripts/hip.fun.R")
pig <- readRDS(file = paste0(inputdir, "OE.liver.all.seurat.filter.final.rds"))
gp_colors <- c("#4472C4", "#70AD47", "#E0BC8D", "#FFC000", "#C83200") %>% 
          setNames(., c("Hepatocytes", "Cholangiocytes", "Endothelial", "Stellate", "immune"))
DimFig(pig, group.by = "group", cols = gp_colors, file_name = "SF1.Liver.lowreso.UMAP", plot.scale = 0.9)



cls_cols <- setNames(c("#6650A2", "#FFA500", "#E1BD8E", "#a31d0e", "#ffcac9", "#9b0104", "#ffc4ca", "#DFD4E9"), 
	c("Hepatocytes", "Stellate", "Cholangiocytes", "myeloid", "B cells", "Plasma cells", "NK_NKT_T", "Endothelial"))
DimFig(pig, group.by = "newcls", cols = cls_cols, file_name = "SF1.Liver.highreso.UMAP", plot.scale = 0.9)








