source("~/Documents/Scripts/hyper.fun.R")

# read in data
data_list <- list(hipo = readRDS(paste0(data, 'hipo.sampled.rds')) %>% 
                    subset(., group %in% "Micro"))

cr_list <- list(infl= read.table(paste0(genelist, 'microglia_inflammatory.csv'), sep = ",", header = F)[,1])
cr_list <- lapply(cr_list, function(z){ z[!is.na(z) & z != ""]})
cls_order <- list(hipo = c("h0", "h1", "h7", "ecmo", "BE"))


# background genes
expg_list <- lapply(data_list, function(seu) {
  exp_genes <-  lapply(levels(as.factor(seu@meta.data$condition)), function(x) {
    subseu <- subset(seu, condition == x)
    genes <- rownames(subseu)[Matrix::rowMeans(subseu$RNA@data !=0) >= 0]
    return(genes)
  })  %>% 
    unlist() %>%
    unique() 
  return(exp_genes)
})

## Extract the cluster markers
mar_list <- lapply(names(data_list), function(ii) {
  mar_fres <- readRDS(paste0(markers, 'Mar_res.subsampled.hipo.micro.rds')) %>%
    subset(pct.1 >= 0.15 & ratio_fc >= 1.1 & avg_log2FC >= 0.15) %>%
    group_by(cluster) %>%
    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
    subset(p_val_adj <= 0.01) 
  allcls <- levels(as.factor(data_list[[ii]]@meta.data$condition))
  mars <- lapply(allcls, function(x) mar_fres$gene[mar_fres$cluster == x]) %>%
    setNames(., allcls)
  return(mars)
}) %>%
  setNames(., names(data_list))


plist <- lapply(names(cr_list), function(sp) {
  pp <- lapply(names(expg_list), function(dataset) {
    exp_genes <- expg_list[[dataset]]
    mars <- mar_list[[dataset]]
    cr_genes <- lapply(cr_list[[sp]], function(x) grep(paste0("^", x, "$"), exp_genes, ignore.case = TRUE, value = TRUE)) %>%
      unlist() %>% unique()
    plot_data <- hp_test(all_genes = exp_genes, marker_list = mars, test_genes = cr_genes) %>%
      data.frame(., stringsAsFactors = FALSE) %>%
      setNames(., "p_val" ) %>%
      rownames_to_column("cluster") %>%
      mutate(mlogp = -log10(p_val))
    plot_data$mlogp[plot_data$mlogp >= 5] <- 5 ## maximum set as 5
    print(min(plot_data$mlogp))
    tt <- paste0(switch(dataset, hipo = "Microglia"), " - ",
                 switch(sp, infl = "pro-inflammatory"))
    p <- ggplot(plot_data, aes_string(x = "cluster", y = "mlogp", fill = "cluster")) +
      geom_bar(size = 0.5, width = 0.8, alpha = 0.9,  stat = "identity") +
      theme(legend.position = "none") +
      theme_bw() +
      labs(title = tt, x = "cluster", y = "-log10(p)") +
      theme_bw() +
      scale_x_discrete(limits = cls_order[[dataset]])+
      scale_y_continuous(limits = c(0, 5)) +
      geom_hline(yintercept = -log10(0.05), size = 0.25, linetype = "dashed")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 7),
            legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
  })
}) %>%
  do.call(c, .)

pdf("Hipo_micro_enrichment.pdf", width = 4, height = 4)
plot_grid(plotlist = plist, nrow = 1, ncol = 1, align = "h")
dev.off()
