source("~/Documents/Scripts/hyper.fun.R")

# read in data
data_list <- list(liver = readRDS(paste0(data, 'liver.sampled.rds')) %>% 
                    subset(., cluster %in% "Hepatocytes"))

cr_list <- list(DNArepair= read.table(paste0(genelist, 'pos_reg_of_dna_repair_GO0045739.txt'), sep = "\t", header = F)[,1],
                Cytoskelet = read.table(paste0(genelist, 'pos_reg_of_cyto_org_GO0051495.txt'), sep = "\t",header = F)[,1],
                Apoptosis= read.table(paste0(genelist, 'negative_reg_of_apop_GO0043066.txt'), sep = "\t", header = F)[,1],
                ATP= read.table(paste0(genelist, 'ATP_met_proc_GO0046034.txt'), sep = "\t", header = F)[,1],
                infl= read.table(paste0(genelist, 'liver_function_markers.csv'), sep = ",", header = F)[,3])

cr_list <- lapply(cr_list, function(z){ z[!is.na(z) & z != ""]})
cls_order <- list(liver = c("h0", "h1", "h7", "ecmo", "BE"))

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
  mar_fres <- readRDS(paste0(markers, 'Mar_res.subsampled.liver.rds')) %>%
    subset(pct.1 >= 0.25 & ratio_fc >= 1.1 & avg_log2FC >= 0.35) %>%
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
    tt <- paste0(switch(dataset, liver = "Hepatocytes"), " - ",
                 switch(sp, 
                        Apoptosis="Neg reg of apoptosis",
                        DNArepair = "Pos. reg. of DNA repair",
                        Cytoskelet = "Positive reorg. of cytoskelet",
                        ATP = "ATP metabolic process",
                        infl="Positive acute phase reactants"))
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
pdf("Liver_enrichment.pdf", width = 10, height = 3)
plot_grid(plotlist = plist, nrow = 1, ncol = 5, align = "h")
dev.off()