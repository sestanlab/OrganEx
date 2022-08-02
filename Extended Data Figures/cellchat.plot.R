
load("../liver.rda")

### Identify and visualize the conserved and context-specific signaling pathways
# overall information flow of each signaling pathway
# comparing OrganEx with other conditions
require(gridExtra)
gg<- list(gg1=rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, comparison = c(5,1), color.use = c("#ffa600", "#003f5c")),
          gg2=rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, comparison = c(5,2),color.use = c("#ffa600", "#003f5c")),
          gg3=rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, comparison = c(5,3),color.use = c("#ffa600", "#003f5c")),
          gg4=rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, comparison = c(5,4),color.use = c("#ffa600", "#79a695")))
pdf("Information_flow.pdf", width = 25, height = 10)
plot_grid(plotlist = gg, ncol=4,nrow = 1,align = "h") %>% print()
dev.off()

### Overall signaling associated with each cell population
#heatmap
library(ComplexHeatmap)
i = 1
pathway.union <- Reduce(union, list(liver.cc.list[[i]]@netP$pathways, liver.cc.list[[i+1]]@netP$pathways, liver.cc.list[[i+2]]@netP$pathways,liver.cc.list[[i+3]]@netP$pathways,liver.cc.list[[i+4]]@netP$pathways))

ht1 = netAnalysis_signalingRole_heatmap(liver.cc.list[[i]], pattern = "all", signaling = pathway.union, title = names(liver.cc.list)[i], width = 8, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(liver.cc.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(liver.cc.list)[i+1], width = 8, height = 20, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(liver.cc.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(liver.cc.list)[i+2], width = 8, height = 20, color.heatmap = "OrRd")
ht4 = netAnalysis_signalingRole_heatmap(liver.cc.list[[i+3]], pattern = "all", signaling = pathway.union, title = names(liver.cc.list)[i+3], width = 8, height = 20, color.heatmap = "OrRd")
ht5 = netAnalysis_signalingRole_heatmap(liver.cc.list[[i+4]], pattern = "all", signaling = pathway.union, title = names(liver.cc.list)[i+4], width = 8, height = 20, color.heatmap = "OrRd")
pdf("Overall_signaling.pdf", width = 20.5, height = 10)
draw(ht1+ht2+ht3+ht4+ht5, ht_gap = unit(0.5, "cm"))
dev.off()
