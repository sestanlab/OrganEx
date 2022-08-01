source("~/project/HIP/scripts/hip.fun.R")
library(Augur)
library(ggpubr)


aug_list1 <- readRDS(file = paste0(inputdir, "Augur.res.combine.rds"))
aug_list2 <- readRDS(file = paste0(inputdir, "Augur.res.combine.sup.rds"))
aug_list <- c(aug_list1, aug_list2)
hip <- readRDS(file = paste0(inputdir, "Augur.HIP.seu.rds"))
hip@meta.data$seurat_clusters <- as.character(hip@meta.data$seurat_clusters)



## For each pair, generate a score
all_ctps <- as.character(0:max(as.numeric(hip@meta.data$seurat_clusters)))
pairs <- c("BE.ecmo", "BE.h0", "BE.h1", "BE.h7", "h0.h1", "h0.h7", "ecmo.h0", "ecmo.h1", "ecmo.h7", "h1.h7")[1:6]


## Bar plots showing major groups
cell2gp <- function(x) {
    yy <- table(x)[table(x)== max(table(x))] %>% names()
    return(yy)
}
cls2gp <- hip@meta.data %>%
        group_by(seurat_clusters) %>%
        summarize(group = cell2gp(group)) %>%
        mutate(seurat_clusters = as.character(seurat_clusters)) %>%
        column_to_rownames("seurat_clusters") %>%
        as.matrix() %>%
        .[, 1]


box_data <- lapply(pairs, function(pair) {
    p1 <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
    p2 <- strsplit(pair, ".", fixed = TRUE)[[1]][2]


    ## Combine the AUC results from donor-pairs
    all_ctps <- as.character(0:23)##aug_list[[pair]][[1]]$AUC$cell_type
    empty <- data.frame(row.names = all_ctps, auc = rep(NA, length(all_ctps)))
    auc_df <- lapply(aug_list[[pair]], function(x) {
        yy <- empty
        yy[x$AUC$cell_type, "auc"] <- x$AUC$auc
        return(yy)
        }) %>%
            do.call(cbind, .) %>%
            as.matrix()

    auc_df[is.na(auc_df)] <- 0.5
    auc_vec <- auc_df %>%
                apply(., 1, function(x) median(x, na.rm = TRUE))

    gp_data <- data.frame(cluster = names(auc_vec),
                            row.names = names(auc_vec), 
                            auc = auc_vec, 
                            group = cls2gp[names(auc_vec)],
                            stringsAsFactors = FALSE,
                            pair = pair)
    return(gp_data)
}) %>%
    do.call(rbind, .)


gp_colors <- c("#6650A2", "#E1BD8E", "#DFD4E9", "#FFA500", "#CD5C5C") %>% 
                setNames(., c("Neuron", "Astro", "OPC Oligo", "Micro", "Vas"))
plist <- lapply(pairs, function(pp) {
    sub_data <- subset(box_data, pair == pp) %>%
                    mutate(group = factor(as.character(group), levels = names(gp_colors)))

    pvals <- sapply(names(gp_colors), function(gp) {
        exp <- sub_data$auc[sub_data$group == gp]
        ref <- sub_data$auc[sub_data$group != gp]
        pval <- wilcox.test(exp, ref, alternative = "greater")$p.val
        pval
        })
	plabel <- sapply(pvals, function(x) {
        if (x <= 0.001){
            return("***")
        } else if (x <= 0.01 & x > 0.001){
            return("**")
        } else if (x <= 0.05 & x > 0.01){
            return("*")
        } else {
            return("")
        }
        })


    anno.data <- data.frame(group = names(pvals), row.names = names(pvals),
                            pval = plabel, 
                            yaxis = 1,
                            stringsAsFactors = FALSE)

    p <- ggplot(sub_data, aes_string(x = "group", y = "auc", fill = "group")) +
          geom_boxplot(size = 0.2, color = "black", outlier.shape = NA) +
          geom_text(data = anno.data, aes_string(x = "group", y = "yaxis", label = "pval"), color = "black", size = 6) +
          scale_fill_manual(values = gp_colors) +
          theme_classic() + 
          labs(title = pp) +
          theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = rel(0.8)), legend.position = "bottom")
    return(p)
    })


pdf(paste0(outputdir, "Augur.HIP.BOX.pdf"), width = 15, height = 10)
plot_grid(plotlist = plist, nrow = 2, ncol = 3) %>% print()
dev.off()






avg_data <- box_data %>%
			mutate(group = factor(as.character(group), levels = names(gp_colors)))
pvals <- sapply(names(gp_colors), function(gp) {
        exp <- avg_data$auc[avg_data$group == gp]
        ref <- avg_data$auc[avg_data$group != gp]
        pval <- wilcox.test(exp, ref, alternative = "greater")$p.val
        pval
        })
plabel <- sapply(pvals, function(x) {
        if (x <= 0.001){
            return("***")
        } else if (x <= 0.01 & x > 0.001){
            return("**")
        } else if (x <= 0.05 & x > 0.01){
            return("*")
        } else {
            return("")
        }
        })

anno.data <- data.frame(group = names(pvals), row.names = names(pvals),
                            pval = plabel, 
                            yaxis = 1,
                            stringsAsFactors = FALSE)
p <- ggplot(avg_data, aes_string(x = "group", y = "auc", fill = "group")) +
          geom_boxplot(size = 0.2, color = "black", outlier.shape = NA) +
          geom_text(data = anno.data, aes_string(x = "group", y = "yaxis", label = "pval"), color = "black", size = 6) +
          scale_fill_manual(values = gp_colors) +
          theme_classic() + 
          theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = rel(0.8)), legend.position = "bottom")
pdf(paste0(outputdir, "Augur.HIP.BOX.summarized.pdf"), width = 4, height = 4)
print(p)
dev.off()




plot_data <- lapply(pairs, function(pair) {
	p1 <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
	p2 <- strsplit(pair, ".", fixed = TRUE)[[1]][2]


	## Combine the AUC results from donor-pairs
	all_ctps <- as.character(0:23)##aug_list[[pair]][[1]]$AUC$cell_type
	empty <- data.frame(row.names = all_ctps, auc = rep(NA, length(all_ctps)))
	auc_df <- lapply(aug_list[[pair]], function(x) {
		yy <- empty
		yy[x$AUC$cell_type, "auc"] <- x$AUC$auc
		return(yy)
		}) %>%
			do.call(cbind, .) %>%
			as.matrix()

	## Use the median AUC as the auc results
	auc_vec <- auc_df %>%
				apply(., 1, function(x) median(x, na.rm = TRUE))

	pair_cells <- colnames(hip)[hip@meta.data$condition %in% c(p1, p2)]
	df <- data.frame(cells = pair_cells,
					group = hip@meta.data[pair_cells, "group"],
					condition = hip@meta.data[pair_cells, "condition"],
					auc = auc_vec[hip@meta.data[pair_cells, "seurat_clusters"]],
					xaxis = hip$umap@cell.embeddings[pair_cells, 1], 
					yaxis = hip$umap@cell.embeddings[pair_cells, 2], 
					pair = pair,
					stringsAsFactors = FALSE)
	return(df)
	})  %>%
		do.call(rbind, .)


plot_data$auc[plot_data$auc < 0.5] <- 0.5
plist <- lapply(pairs, function(pp) {
	p <- ggplot(subset(plot_data, pair == pp), aes_string(x = "xaxis", y = "yaxis", color = "auc")) +
          geom_point(size = 0.01) +
          scale_color_gradientn(colors = c("#dbd8e3", "#6650A2"), limits = c(0.5, 1)) +
          theme_classic() + 
          labs(title = pp) +
          theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "bottom")
    return(p)
	})


jpeg(paste0(outputdir, "Augur.HIP.UMAP.v2.jpeg"), width = 15, height = 11, units = "in", res = 300)
plot_grid(plotlist = plist, nrow = 2, ncol = 3) %>% print()
dev.off() 





