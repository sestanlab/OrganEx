## Calculate the marker res
if (TRUE){
  data_list = data_list
  
  for (ii in names(data_list)){
    seu <- data_list[[ii]]
    all_cls  <- levels(as.factor(seu@meta.data$cluster))
    Idents(seu) <- "cluster"
    marres <- lapply(all_cls, function(cls) {
      aa <- FindMarkers(seu, min.pct = 0.1, ident.1 =cls, ident.2 = NULL, logfc.threshold = 0.25, only.pos = TRUE) %>%
        rownames_to_column("gene") %>%
        mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
        mutate(cluster = cls)
      return(aa)
    }) %>%
      do.call(rbind, .)
    saveRDS(marres, file = paste0(outputdir, "Marker_genes.", ii, ".rds"))
  }
}

markers_list <- list(hipo.mar = readRDS("markers/Marker_genes.hipo.hres.rds"),
                     liver.mar = readRDS("markers/Marker_genes.liver.hres.rds"),
                     kidney.mar = readRDS("markers/Marker_genes.kidney.hres.rds"),
                     heart.mar = readRDS("markers/Marker_genes.heart.hres.rds"))

markers_list.top10 <- lapply(markers_list, function(x) {
  top.10.x <- x %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
})


pdf('heatmaps/Hipo.hres.Heatmap.pdf', height = 15, width = 15)
DoHeatmap(data_list$hipo.hres, group.by = 'cluster', slot = "scale.data", features = markers_list.top10$hipo.mar$gene)

pdf('heatmaps/Liver.hres.Heatmap.pdf', height = 15, width = 15)
DoHeatmap(data_list$liver.hres, group.by = 'cluster', slot = "scale.data", features = markers_list.top10$liver.mar$gene)
dev.off()

pdf('heatmaps/Kidney.hres.Heatmap.pdf', height = 15, width = 15)
DoHeatmap(data_list$kidney.hres, group.by = 'cluster', slot = "scale.data", features = markers_list.top10$kidney.mar$gene)
dev.off()

pdf('heatmaps/Heart.hres.Heatmap.pdf', height = 15, width = 15)
DoHeatmap(data_list$heart.hres, group.by = 'cluster', slot = "scale.data", features = markers_list.top10$heart.mar$gene) 
dev.off()
