## Calculate the marker res
if (TRUE){
  heart.all <- readRDS(paste0(inputdir, "heart.sampled.rds")) %>%
    subset(., cluster %in% 'Cardiomyocytes')
  heart.functional <- readRDS(paste0(inputdir, "heart.sampled.rds")) %>%
    subset(., condition %in% c('h0','ecmo','BE') & cluster %in% 'Cardiomyocytes')
  hipo.neurons <- readRDS(paste0(inputdir, "hipo.sampled.rds")) %>%
    subset(., group %in% 'Neuron')
  hipo.astro <- readRDS(paste0(inputdir, "hipo.sampled.rds")) %>%
    subset(., group %in% 'Astro')
  hipo.micro <- readRDS(paste0(inputdir, "hipo.sampled.rds")) %>%
    subset(., group %in% 'Micro')
  liver <- readRDS(paste0(inputdir, "liver.sampled.rds")) %>%
    subset(., cluster %in% 'Hepatocytes')
  kidney <- readRDS(paste0(inputdir, "kidney.sampled.rds")) %>%
    subset(., cluster %in% 'Proximal Tubule')
  
  data_list <- list(heart.all = heart.all, heart.functional = heart.functional, hipo.neurons = hipo.neurons, hipo.astro = hipo.astro, hipo.micro = hipo.micro, liver = liver, kidney = kidney)
  
  for (ii in names(data_list)){
    seu <- data_list[[ii]]
    all_cls  <- levels(as.factor(seu@meta.data$condition))
    Idents(seu) <- "condition"
    marres <- lapply(all_cls, function(cls) {
      aa <- FindMarkers(seu, min.pct = 0.1, ident.1 =cls, ident.2 = NULL, logfc.threshold = 0.1, only.pos = TRUE) %>%
        rownames_to_column("gene") %>%
        mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
        mutate(cluster = cls)
      return(aa)
    }) %>%
      do.call(rbind, .)
    saveRDS(marres, file = paste0("Mar_res.subsampled.", ii, ".rds"))
  }
}


## Calculate the marker res
if (TRUE){
  liver <- readRDS(paste0(inputdir, "liver.sampled.rds")) %>%
    subset(., cluster %in% 'Hepatocytes' & condition %in% c('h0', 'ecmo', 'BE'))
  kidney <- readRDS(paste0(inputdir, "kidney.sampled.rds")) %>%
    subset(., cluster %in% 'Proximal Tubule' & condition %in% c('h0', 'ecmo', 'BE'))
  
  data_list <- list(liver = liver, kidney = kidney)
  
  for (ii in names(data_list)){
    seu <- data_list[[ii]]
    all_cls  <- levels(as.factor(seu@meta.data$condition))
    Idents(seu) <- "condition"
    marres <- lapply(all_cls, function(cls) {
      aa <- FindMarkers(seu, min.pct = 0.1, ident.1 =cls, ident.2 = NULL, logfc.threshold = 0.1, only.pos = TRUE) %>%
        rownames_to_column("gene") %>%
        mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
        mutate(cluster = cls)
      return(aa)
    }) %>%
      do.call(rbind, .)
    saveRDS(marres, file = paste0("Mar_res.subsampled.", ii, ".rds"))
  }
}