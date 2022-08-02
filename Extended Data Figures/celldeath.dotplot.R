data_list <- list(hipo = readRDS("data/OE.HIP.all.seurat.filter.final.rds") %>% 
                    subset(., group %in% "Neuron"),
                  liver = readRDS("data/OE.liver.all.seurat.filter.final.rds") %>% 
                    subset(., cluster %in% "Hepatocytes"),
                  kidney = readRDS("data/OE.kidney.all.seurat.filter.final.rds") %>% 
                    subset(., cluster %in% "Proximal Tubule"),
                  heart = readRDS("data/OE.heart.all.seurat.filter.final.rds") %>% 
                    subset(., cluster %in% "Cardiomyocytes"))

gene.list <- read.csv('NEW_ANALYSIS/go_terms/cell_death_pathways.csv', header = TRUE)
gene.list <- lapply(gene.list, function(z){ z[!is.na(z) & z != ""]})
all.genes <- unique(unlist(gene.list))

b<-'condition'
List<-list()
for (i in 1:length(data_list)) {
  for (ii in 1:length(b)) {
    print(i);print(ii)
    plot=DotPlot(object = data_list[[i]],group.by=b[ii], features = all.genes, cols = c("#90AFC5", "#A43820")) + RotatedAxis() + scale_size(range = c(5,15)) + ggtitle(paste0("Death pathways - ", names(data_list)[i]))
    List[[length(List)+1]]<-plot 
    for(z in 1:length(List)){
      List[[z]]$data$features.plot <- factor(List[[z]]$data$features.plot, levels = c(gene.list$intrinsic, gene.list$p53, gene.list$extrinsic, gene.list$pyro,gene.list$ferro, gene.list$necro))
      List[[z]]$data$id <- factor(List[[z]]$data$id, levels= c('BE', 'ecmo', 'h7', 'h1', 'h0'))
    }
  }
}

pdf(paste0(outputdir, "DotPlots.pdf"), width = 40, height = 40)
plot_grid(plot_grid(plotlist = List, ncol = 1, nrow = 4)) %>% print()
dev.off()
