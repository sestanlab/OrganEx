## Use the Augur to find the most different cluster among conditions
args <- commandArgs(trailingOnly = TRUE)
source("~/project/HIP/scripts/hip.fun.R")
library(Augur)
library(ggpubr)


seu_name <- "final"
file_name <- "Seurat.withp46"


pairs <- c("BE.ecmo", "BE.h0", "BE.h1", "BE.h7", "ecmo.h0", "ecmo.h1", "ecmo.h7", "h0.h1", "h0.h7", "h1.h7")
pair <- pairs[as.numeric(args[1])]



if (FALSE){
	## Find Clusters
	hip <- readRDS(file = paste0(inputdir, "OE.HIP.all.seurat.filter.", "final", ".rds"))
	hip <- hip[, hip@meta.data$samplename != "p46_BE"]
	hip <- FindNeighbors(hip, dims = 1:30, k.param = 25) %>%
	                        FindClusters(., resolution = 0.5, n.iter = 20, random.seed = 42)
	DimFig(hip, group.by = "seurat_clusters", plot.scale = 0.9, file_name = "Augur.HIP")

	##Rationale to merge clusters
	##avg <- log(AverageExpression(hip, features = rownames(hip$pca@feature.loadings))$RNA + 1)
	##cor(avg[, c("0", "1", "26", "27", "28", "29", "30", "31")], method = "p")
	hip@meta.data$seurat_clusters[hip@meta.data$seurat_clusters %in% c("26", "27", "28", "29", "30", "31")] <- "0"
	hip@meta.data$seurat_clusters[hip@meta.data$seurat_clusters %in% c("24")] <- "14"
	hip@meta.data$seurat_clusters[hip@meta.data$seurat_clusters %in% c("25")] <- "7"
	DimFig(hip, group.by = "seurat_clusters", plot.scale = 0.9, file_name = "Augur.HIP.merge")
	saveRDS(hip, file = paste0(inputdir, "Augur.HIP.seu.rds"))
}



## Calculate the Augur results
p1 <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
p2 <- strsplit(pair, ".", fixed = TRUE)[[1]][2]

hip <- readRDS(file = paste0(inputdir, "Augur.HIP.seu.rds"))
hip@meta.data$seurat_clusters <- as.character(hip@meta.data$seurat_clusters)


sp1 <- hip@meta.data$samplename[hip@meta.data$condition == p1] %>% unique()
sp2 <- hip@meta.data$samplename[hip@meta.data$condition == p2] %>% unique()
subseu <- hip[, hip@meta.data$condition %in% c(p1, p2)]
subseu <- subseu[rownames(hip$pca@feature.loadings), ]


all_aug <- lapply(sp1, function(xx){
	aug_cbn <- lapply(sp2, function(yy) {
		sp_seu <- subseu[, subseu@meta.data$samplename %in% c(xx, yy)]
		sub_aug <- calculate_auc(sp_seu$RNA@data, meta = sp_seu@meta.data, cell_type_col = "seurat_clusters", label_col = "condition", n_threads = 8)
		return(sub_aug)
		}) %>%
			setNames(., paste0(xx, ".", sp2))
	return(aug_cbn)
	}) %>% 
		do.call(c, .)
saveRDS(all_aug, file = paste0(inputdir, "Augur.res.", pair, ".rds"))




## Combine the Augur results
if (FALSE){
	pairs <- c("BE.ecmo", "BE.h0", "BE.h1", "BE.h7", "ecmo.h0", "ecmo.h1", "ecmo.h7", "h0.h1", "h0.h7", "h1.h7")

	aug_list <- lapply(pairs, function(pair) {
		readRDS(file = paste0(inputdir, "Augur.res.", pair, ".rds"))
		}) %>%
			setNames(., pairs)
	saveRDS(aug_list, file = paste0(inputdir, "Augur.res.combine.rds"))
}



