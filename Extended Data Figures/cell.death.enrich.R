# calculate markers
data_list <- list(hipo = readRDS(paste0(inputdir, "OE.HIP.all.seurat.filter.final.rds")) %>% 
                    subset(., group %in% "Neuron" & condition %in% c('h0','ecmo', 'BE')),
                  liver = readRDS(paste0(inputdir, "OE.liver.all.seurat.filter.final.rds")) %>% 
                    subset(., cluster %in% "Hepatocytes" & condition %in% c('h0','ecmo', 'BE')),
                  kidney = readRDS(paste0(inputdir, "OE.kidney.all.seurat.filter.final.rds")) %>% 
                    subset(., cluster %in% "Proximal Tubule" & condition %in% c('h0','ecmo', 'BE')),
                  heart = readRDS(paste0(inputdir, "OE.heart.all.seurat.filter.final.rds")) %>% 
                    subset(., cluster %in% "Cardiomyocytes" & condition %in% c('h0','ecmo', 'BE')))

for (ii in names(data_list)){
  seu <- data_list[[ii]]
  all_cls  <- levels(as.factor(seu@meta.data$condition))
  Idents(seu) <- 'condition'
  marres <- lapply(all_cls, function(cls) {
    aa <- FindMarkers(seu, min.pct = -1, ident.1 =cls, ident.2 = NULL, logfc.threshold = -1) %>%
      rownames_to_column("gene") %>%
      mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
      mutate(cluster = cls)
    return(aa)
  }) %>%
    do.call(rbind, .)
  saveRDS(marres, file = paste0("Mar_res.", ii, ".rds"))
}

#### hipo
hipo <- readRDS("../../markers/Mar_res.hipo.rds")
hipo <- split(hipo, hipo$cluster)
hipo <- lapply(hipo, function(x) mutate(x, p_val_adj = p.adjust(p_val, method = "fdr")) %>% mutate(mlogp = -log10(p_val_adj)) %>% group_by(cluster) %>% subset(gene %in% all.genes))

for(z in 1:length(hipo)){
  hipo[[z]]$mlogp[hipo[[z]]$mlogp == "Inf" & hipo[[z]]$avg_log2FC > 0]  <- 100
  hipo[[z]]$mlogp[hipo[[z]]$avg_log2FC < 0]  <- 0
  hipo[[z]]$mlogp[hipo[[z]]$mlogp >= 100 & hipo[[z]]$avg_log2FC > 0]  <- 100
}

for (i in 1:length(hipo)){
  df_i <- hipo[[i]]
  pp <-  
    df_i %>%
    ggplot(aes(x = reorder(gene, -mlogp), y = mlogp)) + 
    geom_bar(size = 0.5, width = 0.8, alpha = 0.9,  stat = "identity", fill = "#003f5c") +
    theme(legend.position = "none") + 
    theme_bw() + 
    scale_y_continuous(limits = c(0, 100)) +
    geom_hline(yintercept = -log10(0.01), size = 0.25, linetype = "dashed")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 7), legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  ggsave(pp, file=paste0(names(hipo)[i],"_barplot.pdf"))
}

### heart
heart <- readRDS("../../markers/Mar_res.heart.rds")
heart <- split(heart, heart$cluster)
heart <- lapply(heart, function(x) mutate(x, p_val_adj = p.adjust(p_val, method = "fdr")) %>% mutate(mlogp = -log10(p_val_adj)) %>% group_by(cluster) %>% subset(gene %in% all.genes))

for(z in 1:length(heart)){
  heart[[z]]$mlogp[heart[[z]]$mlogp == "Inf" & heart[[z]]$avg_log2FC > 0]  <- 100
  heart[[z]]$mlogp[heart[[z]]$avg_log2FC < 0]  <- 0
  heart[[z]]$mlogp[heart[[z]]$mlogp >= 100 & heart[[z]]$avg_log2FC > 0]  <- 100
}

for (i in 1:length(heart)){
  df_i <- heart[[i]]
  pp <-  
    df_i %>%
    ggplot(aes(x = reorder(gene, -mlogp), y = mlogp)) + 
    geom_bar(size = 0.5, width = 0.8, alpha = 0.9,  stat = "identity", fill = "#003f5c") +
    theme(legend.position = "none") + 
    theme_bw() + 
    scale_y_continuous(limits = c(0, 100)) +
    geom_hline(yintercept = -log10(0.01), size = 0.25, linetype = "dashed")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 7), legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  ggsave(pp, file=paste0(names(heart)[i],"_barplot.pdf"))
}

### liver hepatocytes
liver <- readRDS("../../markers/Mar_res.liver.rds")
liver <- split(liver, liver$cluster)
liver <- lapply(liver, function(x) mutate(x, p_val_adj = p.adjust(p_val, method = "fdr")) %>% mutate(mlogp = -log10(p_val_adj)) %>% group_by(cluster) %>% subset(gene %in% all.genes))

for(z in 1:length(liver)){
  liver[[z]]$mlogp[liver[[z]]$mlogp == "Inf" & liver[[z]]$avg_log2FC > 0]  <- 100
  liver[[z]]$mlogp[liver[[z]]$avg_log2FC < 0]  <- 0
  liver[[z]]$mlogp[liver[[z]]$mlogp >= 100 & liver[[z]]$avg_log2FC > 0]  <- 100
}

for (i in 1:length(liver)){
  df_i <- liver[[i]]
  pp <-  
    df_i %>%
    ggplot(aes(x = reorder(gene, -mlogp), y = mlogp)) + 
    geom_bar(size = 0.5, width = 0.8, alpha = 0.9,  stat = "identity", fill = "#003f5c") +
    theme(legend.position = "none") + 
    theme_bw() + 
    scale_y_continuous(limits = c(0, 100)) +
    geom_hline(yintercept = -log10(0.01), size = 0.25, linetype = "dashed")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 7), legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  ggsave(pp, file=paste0(names(liver)[i],"_barplot.pdf"))
}

### kidney
kidney <- readRDS("../../markers/Mar_res.kidney.rds")
kidney <- split(kidney, kidney$cluster)
kidney <- lapply(kidney, function(x) mutate(x, p_val_adj = p.adjust(p_val, method = "fdr")) %>% mutate(mlogp = -log10(p_val_adj)) %>% group_by(cluster) %>% subset(gene %in% all.genes))

for(z in 1:length(kidney)){
  kidney[[z]]$mlogp[kidney[[z]]$mlogp == "Inf" & kidney[[z]]$avg_log2FC > 0]  <- 100
  kidney[[z]]$mlogp[kidney[[z]]$avg_log2FC < 0]  <- 0
  kidney[[z]]$mlogp[kidney[[z]]$mlogp >= 100 & kidney[[z]]$avg_log2FC > 0]  <- 100
}

for (i in 1:length(kidney)){
  df_i <- kidney[[i]]
  pp <-  
    df_i %>%
    ggplot(aes(x = reorder(gene, -mlogp), y = mlogp)) + 
    geom_bar(size = 0.5, width = 0.8, alpha = 0.9,  stat = "identity", fill = "#003f5c") +
    theme(legend.position = "none") + 
    theme_bw() + 
    scale_y_continuous(limits = c(0, 100)) +
    geom_hline(yintercept = -log10(0.01), size = 0.25, linetype = "dashed")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 7), legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  ggsave(pp, file=paste0(names(kidney)[i],"_barplot.pdf"))
}
