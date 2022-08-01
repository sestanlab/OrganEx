## Use the Augur to find the most different clusters among conditions
args <- commandArgs(trailingOnly = TRUE)
source("~/project/HIP/scripts/hip.fun.R")
library(Augur)
library(ggpubr)


file_name <- "Seurat"
pairs <- c("BE.ecmo", "BE.h0", "BE.h1", "BE.h7", "ecmo.h0", "ecmo.h1", "ecmo.h7", "h0.h1", "h0.h7", "h1.h7")
pair <- pairs[as.numeric(args[1])]



## Calculate the Augur results
p1 <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
p2 <- strsplit(pair, ".", fixed = TRUE)[[1]][2]

hip <- readRDS(file = paste0(inputdir, "OE.heart.all.seurat.filter.final.rds"))
hip@meta.data$seurat_clusters <- as.character(hip@meta.data$seurat_clusters)


sp1 <- hip@meta.data$samplename[hip@meta.data$condition == p1] %>% unique()
sp2 <- hip@meta.data$samplename[hip@meta.data$condition == p2] %>% unique()
subseu <- hip[, hip@meta.data$condition %in% c(p1, p2)]

exp_genes <- rownames(subseu)[Matrix::rowSums(subseu$RNA@data != 0) >= 30]
print(length(exp_genes))
subseu <- subseu[exp_genes, ]


all_aug <- lapply(sp1, function(xx){
	aug_cbn <- lapply(sp2, function(yy) {
		sp_seu <- subseu[, subseu@meta.data$samplename %in% c(xx, yy)]
		sub_aug <- calculate_auc(sp_seu$RNA@data, meta = sp_seu@meta.data, cell_type_col = "seurat_clusters", label_col = "condition", n_threads = 8)
		print(paste0("Finish sample pair: ", xx, " - ", yy))
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



