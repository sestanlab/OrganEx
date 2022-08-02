## Calculate the marker res
if (TRUE){
  heart <- readRDS(paste0(inputdir, "OE.heart.all.seurat.filter.final.rds")) %>%
    subset(., cluster %in% 'Cardiomyocytes' & condition %in% c('h0', 'h1', 'h7'))
    hipo <- readRDS(paste0(inputdir,"OE.HIP.all.seurat.filter.final.rds")) %>%
    subset(., group %in% 'Neuron' & condition %in% c('h0', 'h1', 'h7'))
    liver <- readRDS(paste0(inputdir,"OE.liver.all.seurat.filter.final.rds")) %>%
    subset(., cluster %in% 'Hepatocytes' & condition %in% c('h0', 'h1', 'h7'))
    kidney <- readRDS(paste0(inputdir,"OE.kidney.all.seurat.filter.final.rds")) %>%
    subset(., cluster %in% 'Proximal Tubule' & condition %in% c('h0', 'h1', 'h7'))

  data_list <- list(hipo=hipo, heart=heart, liver=liver, kidney=kidney)

  for (ii in names(data_list)){
    seu <- data_list[[ii]]
    all_cls  <- c("h1","h7")
    Idents(seu) <- "condition"
    marres <- lapply(all_cls, function(cls) {
      aa <- FindMarkers(seu, min.pct = 0.1, ident.1 ="h0", ident.2 = cls, logfc.threshold = 0.1, only.pos = FALSE, max.cells.per.ident = 1000) %>%
        rownames_to_column("gene") %>%
        mutate(ratio_fc1 = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
        mutate(ratio_fc2 = (pct.2 + 0.01)/(pct.1 + 0.01)) %>%
        mutate(cluster = cls)
      return(aa)
    }) %>%
      do.call(rbind, .)
    saveRDS(marres, file = paste0("DEG.OEx.", ii, ".rds"))
  }
}
DEG <- list(hipo = readRDS(paste0(genes, "DEG.OEx.hipo.rds")),
            heart = readRDS(paste0(genes, "DEG.OEx.heart.rds")),
            liver = readRDS(paste0(genes, "DEG.OEx.liver.rds")),
            kidney = readRDS(paste0(genes, "DEG.OEx.kidney.rds")))

DEG.pos <- lapply(DEG, function(x)filter(x, avg_log2FC >= 0.1 & ratio_fc1 > 1) %>%
                   mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>% subset(p_val_adj < 0.05) %>%
                    group_by(cluster) %>%
                    top_n(800, wt = avg_log2FC))
DEG.neg <- lapply(DEG, function(x)filter(x, avg_log2FC <= -0.1 & ratio_fc2 > 1) %>%
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>% subset(p_val_adj < 0.05) %>%
                    group_by(cluster) %>%
                    top_n(-800, wt = avg_log2FC))

DEG <- list(hipo.h1.genes=hipo.h1.genes, hipo.h7.genes=hipo.h7.genes,
			heart.h1.genes=heart.h1.genes, heart.h7.genes=heart.h7.genes,
			kidney.h1.genes=kidney.h1.genes, kidney.h7.genes=kidney.h7.genes,
			liver.h1.genes=liver.h1.genes, liver.h7.genes=liver.h7.genes)

GOres = list() 
length(GOres) = length(DEG) 
names(GOres) = names(DEG) # name elements the same as in DEGs list
for(i in 1:length(DEG)) { # loop over the list of DEGs
	GOres[[i]] = enrichGO(gene = DEG[[i]],
		OrgDb = org.Ss.eg.db, # Sus Scrofa
		ont = "BP", 
		pAdjustMethod = "fdr",
		keyType = 'SYMBOL',
		pvalueCutoff = 0.05,
		qvalueCutoff = 0.2,
		maxGSSize = 1000,
		readable = F)
}

results <- list()
length(results) <- length(GOres)
names(results) = names(GOres)

for(i in 1:length(GOres)) { 
	results[[i]] <- GOres[[i]]@result
}
results <- lapply(results, function(x) mutate(x, mlogp = -log10(p.adjust)) %>%
			top_n(-15, wt = p.adjust))
  
# plot
for (i in 1:length(results)){
	results[[i]]$mlogp[results[[i]]$mlogp >= 5] <- 5
	df_i <- results[[i]]
	pp <-  
			df_i %>%
			ggplot(aes(x = mlogp, y = reorder(Description, -mlogp))) + 
			geom_bar(width = 0.5, stat = "identity", fill = "#A7A9AC") +
			theme(legend.position = "none") +
			scale_x_continuous(limits = c(0, 5)) + 
			geom_vline(xintercept = -log10(0.05), size = 0.25, linetype = "dashed") +
			theme_classic() + ggtitle(names(results)[i]) +
			theme(legend.position = "none",
			          title = element_text(size = 3),
			          axis.title = element_blank(),
			          axis.text.x = element_text(size = 3),
			          axis.ticks.length.x = unit(0.7, "mm"),
			          axis.ticks.x = element_line(size = 0.4),
			          axis.text.y = element_text(size = 5),
			          axis.ticks.y = element_blank(), axis.line.y = element_blank())
ggsave(pp, file=paste0(output.dir, names(results)[i],"_GO_barplot.pdf"), width = 8, height = 3)
}



