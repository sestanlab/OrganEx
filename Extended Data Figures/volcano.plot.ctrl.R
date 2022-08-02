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
      aa <- FindMarkers(seu, min.pct = -1, ident.1 ="h0", ident.2 = cls, logfc.threshold = -1, only.pos = FALSE, max.cells.per.ident = 1000) %>%
        rownames_to_column("gene") %>%
        mutate(ratio_fc1 = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
        mutate(ratio_fc2 = (pct.2 + 0.01)/(pct.1 + 0.01)) %>%
        mutate(cluster = cls)
      return(aa)
    }) %>%
      do.call(rbind, .)
    saveRDS(marres, file = paste0("DEG.OEx.volcano.", ii, ".rds"))
  }
}

DEG <- list(hipo = readRDS(paste0(genes, "DEG.OEx.volcano.hipo.rds")),
            heart = readRDS(paste0(genes, "DEG.OEx.volcano.heart.rds")),
            liver = readRDS(paste0(genes, "DEG.OEx.volcano.liver.rds")),
            kidney = readRDS(paste0(genes, "DEG.OEx.volcano.kidney.rds")))

DEG.pos <- lapply(DEG, function(x)filter(x, avg_log2FC >= 0.25 & ratio_fc1 > 1) %>%
                    filter(!grepl("rna", gene)) %>% filter(!grepl("LOC", gene)) %>%
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>% subset(p_val_adj < 0.05) %>%
                    group_by(cluster) %>%
                    top_n(15, wt = avg_log2FC))
DEG.neg <- lapply(DEG, function(x)filter(x, avg_log2FC <= -0.25 & ratio_fc2 > 1) %>%
                    filter(!grepl("rna", gene)) %>% filter(!grepl("LOC", gene)) %>%
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>% subset(p_val_adj < 0.05) %>%
                    group_by(cluster) %>%
                    top_n(-15, wt = avg_log2FC))

DEG <- lapply(DEG, function(x) mutate(x, mlogp = -log10(p_val_adj)))
for(z in 1:length(DEG)){
  DEG[[z]]$p_val_adj[DEG[[z]]$p_val_adj <= 6.109245e-51]  <- 6.109245e-51
}


####### hipo
hipo <- split(DEG$hipo, DEG$hipo$cluster)
hipo.deg.pos <- split(DEG.pos$hipo, DEG.pos$hipo$cluster)
hipo.deg.neg <- split(DEG.neg$hipo, DEG.neg$hipo$cluster)

results = list() 
length(results) = length(hipo) 
names(results) = names(hipo) 
# plotting 
for(i in 1:length(hipo)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(hipo[[i]],
                                 lab = hipo[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 selectLab = c(DEG.pos[[i]]$gene, DEG.neg[[i]]$gene),
                                 labSize = 1,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 subtitle = NA,
                                 caption = NA,
                                 xlim = c(-1.5, 1.5),
                                 title = names(hipo[i])) + theme_classic() +
                                 theme(axis.line = element_line(size = 0.5), legend.position = "none",
                                 axis.text = element_text(size = 3), )
}

#export as pdf
for(i in 1:length(results)) { # loop over the results
  pdf(file = paste0(output.dir, "Volcano.frame.hip.h0.vs.", names(results[i]), ".pdf"), height = 8, width = 8)
  plot(results[[i]])
  dev.off()
}

# heart
heart <- split(DEG$heart, DEG$heart$cluster)
heart.deg.pos <- split(DEG.pos$heart, DEG.pos$heart$cluster)
heart.deg.neg <- split(DEG.neg$heart, DEG.neg$heart$cluster)

results = list() 
length(results) = length(heart) 
names(results) = names(heart)
# plotting
for(i in 1:length(heart)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(heart[[i]],
                                 lab = heart[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 selectLab = c(DEG.pos[[i]]$gene, DEG.neg[[i]]$gene),
                                 labSize = 1,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 subtitle = NA,
                                 caption = NA,
                                 xlim = c(-2.1, 2),
                                 title = names(heart[i])) + theme_classic() +
    theme(axis.line = element_line(size = 0.5), legend.position = "none",
          axis.text = element_text(size = 3), )
}

#export as pdf
for(i in 1:length(results)) { # loop over the results
  pdf(file = paste0(output.dir, "Volcano.frame.heart.h0.vs.", names(results[i]), ".pdf"), height = 8, width = 8)
  plot(results[[i]])
  dev.off()
}


# kidney
kidney <- split(DEG$kidney, DEG$kidney$cluster)
kidney.deg.pos <- split(DEG.pos$kidney, DEG.pos$kidney$cluster)
kidney.deg.neg <- split(DEG.neg$kidney, DEG.neg$kidney$cluster)

results = list()
length(results) = length(kidney) 
names(results) = names(kidney)
# plotting axes
for(i in 1:length(kidney)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(kidney[[i]],
                                 lab = kidney[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 selectLab = c(DEG.pos[[i]]$gene, DEG.neg[[i]]$gene),
                                 labSize = 1,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 subtitle = NA,
                                 caption = NA,
                                 xlim = c(-2.5, 2.5),
                                 title = names(kidney[i])) + theme_classic() +
    theme(axis.line = element_line(size = 0.5), legend.position = "none",
          axis.text = element_text(size = 3), )
}

#export as pdf
for(i in 1:length(results)) { # loop over the results
  pdf(file = paste0(output.dir, "Volcano.frame.kidney.h0.vs.", names(results[i]), ".pdf"), height = 8, width = 8)
  plot(results[[i]])
  dev.off()
}


# liver
liver <- split(DEG$liver, DEG$liver$cluster)
liver.deg.pos <- split(DEG.pos$liver, DEG.pos$liver$cluster)
liver.deg.neg <- split(DEG.neg$liver, DEG.neg$liver$cluster)

results = list()
length(results) = length(liver) 
names(results) = names(liver) 
# plotting axes
for(i in 1:length(liver)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(liver[[i]],
                                 lab = liver[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 selectLab = c(DEG.pos[[i]]$gene, DEG.neg[[i]]$gene),
                                 labSize = 1,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 subtitle = NA,
                                 caption = NA,
                                 xlim = c(-4.1, 4),
                                 title = names(liver[i])) + theme_classic() +
    theme(axis.line = element_line(size = 0.5), legend.position = "none",
          axis.text = element_text(size = 3), )
}

#export as pdf
for(i in 1:length(results)) { # loop over the results
  pdf(file = paste0(output.dir, "Volcano.frame.liver.h0.vs.", names(results[i]), ".pdf"), height = 8, width = 8)
  plot(results[[i]])
  dev.off()
}


