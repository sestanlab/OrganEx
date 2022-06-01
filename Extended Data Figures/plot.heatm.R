load('subsampled/OE.seurat.objects.hres.subsampled.RData')

markers_list <- list(hipo.mar = readRDS("markers/Marker_genes.hipo.hres.rds"),
                     liver.mar = readRDS("markers/Marker_genes.liver.hres.rds"),
                     kidney.mar = readRDS("markers/Marker_genes.kidney.hres.rds"),
                     heart.mar = readRDS("markers/Marker_genes.heart.hres.rds"))

markers_list.top10 <- lapply(markers_list, function(x) {
  top.10.x <- x %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
})


pdf('heatmaps/Hipo.hres.Heatmap.pdf', height = 15, width = 15)
DoHeatmap(data_list$hipo.hres, group.by = 'cluster', slot = "scale.data", features = markers_list.top10$hipo.mar$gene)
dev.off()
