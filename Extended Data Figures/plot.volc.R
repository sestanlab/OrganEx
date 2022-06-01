
#read in data
DEG <- list(hipo = readRDS(paste0(genes, "DEG.OEx.volcano.hipo.rds")),
            heart = readRDS(paste0(genes, "DEG.OEx.volcano.heart.rds")),
            liver = readRDS(paste0(genes, "DEG.OEx.volcano.liver.rds")),
            kidney = readRDS(paste0(genes, "DEG.OEx.volcano.kidney.rds")))

for(z in 1:length(DEG)){
  DEG[[z]]$p_val_adj[DEG[[z]]$p_val_adj <= 1e-50]  <- 1e-50
}
DEG <- lapply(DEG, function(x) mutate(x, mlogp = -log10(p_val_adj)))

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



####### hipo
hipo <- split(DEG$hipo, DEG$hipo$cluster)
hipo.deg.pos <- split(DEG.pos$hipo, DEG.pos$hipo$cluster)
hipo.deg.neg <- split(DEG.neg$hipo, DEG.neg$hipo$cluster)

results = list() 
length(results) = length(hipo) 
names(results) = names(hipo) 

for(i in 1:length(hipo)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(hipo[[i]],
                                 lab = hipo[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 labCol = NA,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 cutoffLineCol = NA,
                                 xlim = c(-2, 2),
                                 ylim = c(0, max(-log10(hipo[[i]]$p_val_adj), na.rm = TRUE) + 5),
                                 title = names(hipo[i])) + theme_classic()
}
for(i in 1:length(results)) { # loop over the  results and save jpeg
  jpeg(file = paste0(output.dir, "Volcano.hip.h0.vs.", names(results[i]), ".jpeg"), res=800, width=1400, height=1100, pointsize=10)
  plot(results[[i]])
  dev.off()
}

###### heart
heart <- split(DEG$heart, DEG$heart$cluster)
heart.deg.pos <- split(DEG.pos$heart, DEG.pos$heart$cluster)
heart.deg.neg <- split(DEG.neg$heart, DEG.neg$heart$cluster)

results = list() 
length(results) = length(heart) 
names(results) = names(heart)

for(i in 1:length(heart)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(heart[[i]],
                                 lab = heart[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 labCol = NA,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 cutoffLineCol = NA,
                                 xlim = c(-3.1, 3.1),
                                 ylim = c(0, max(-log10(heart[[i]]$p_val_adj), na.rm = TRUE) + 5),
                                 title = names(heart[i])) + theme_classic()
}

for(i in 1:length(results)) { # loop over the  results and save jpeg
  jpeg(file = paste0(output.dir, "Volcano.heart.h0.vs.", names(results[i]), ".jpeg"), res=800, width=1400, height=1100, pointsize=10)
  plot(results[[i]])
  dev.off()
}

######## kidney
kidney <- split(DEG$kidney, DEG$kidney$cluster)
kidney.deg.pos <- split(DEG.pos$kidney, DEG.pos$kidney$cluster)
kidney.deg.neg <- split(DEG.neg$kidney, DEG.neg$kidney$cluster)

results = list() 
length(results) = length(kidney) 
names(results) = names(kidney) 

for(i in 1:length(kidney)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(kidney[[i]],
                                 lab = kidney[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 labCol = NA,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 cutoffLineCol = NA,
                                 xlim = c(-4.6, 4.6),
                                 ylim = c(0, max(-log10(kidney[[i]]$p_val_adj), na.rm = TRUE) + 5),
                                 title = names(kidney[i])) + theme_classic()
}
for(i in 1:length(results)) { # loop over the  results and save jpeg
  jpeg(file = paste0(output.dir, "Volcano.kidney.h0.vs.", names(results[i]), ".jpeg"), res=800, width=1400, height=1100, pointsize=10)
  plot(results[[i]])
  dev.off()
}


####### liver
liver <- split(DEG$liver, DEG$liver$cluster)
liver.deg.pos <- split(DEG.pos$liver, DEG.pos$liver$cluster)
liver.deg.neg <- split(DEG.neg$liver, DEG.neg$liver$cluster)

results = list() 
length(results) = length(liver) 
names(results) = names(liver)

for(i in 1:length(liver)) { # loop over the list of DEGs
  results[[i]] = EnhancedVolcano(liver[[i]],
                                 lab = liver[[i]]$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 labCol = NA,
                                 col = c("grey81", "grey81", "#E1BD8E", "#CD5C5C"),
                                 pointSize = 0.1,
                                 pCutoff = 0.01,
                                 FCcutoff = 0.25,
                                 cutoffLineCol = NA,
                                 xlim = c(-3.6, 3.6),
                                 ylim = c(0, max(-log10(liver[[i]]$p_val_adj), na.rm = TRUE) + 5),
                                 title = names(liver[i])) + theme_classic()
}
for(i in 1:length(results)) { # loop over the  results and save jpeg
  jpeg(file = paste0(output.dir, "Volcano.liver.h0.vs.", names(results[i]), ".jpeg"), res=800, width=1400, height=1100, pointsize=10)
  plot(results[[i]])
  dev.off()
}