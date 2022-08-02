# heart
heart <- readRDS(paste0(inputdir,'OE.heart.all.seurat.filter.final.rds'))
heart@meta.data$cdt_gp <- paste0(as.character(heart@meta.data$condition), "_", as.character(heart@meta.data$group))
heart.list <- SplitObject(heart, split.by = "cdt_gp")

set.seed(0)
heart.list.sampled <- lapply(X = heart.list, N = 1000, FUN = function(x, N) {
  if(ncol(x) > N){
    x <- x[,sample(ncol(x), N, replace = FALSE)] 
  } else {
    x
  }
  return(x)
})

heart.sampled <- Reduce(merge,heart.list.sampled)
saveRDS(heart.sampled, "heart.sampled.rds")

heart <- readRDS(paste0(inputdir, 'heart.sampled.rds'))
heart.list <- list(h0= createCellChat(object = subset(heart, subset = condition == "h0"), group.by = "group"),
                  h1= createCellChat(object = subset(heart, subset = condition == "h1"), group.by = "group"),
                  h7= createCellChat(object = subset(heart, subset = condition == "h7"), group.by = "group"),
                  ecmo= createCellChat(object = subset(heart, subset = condition == "ecmo"), group.by = "group"),
                  be=createCellChat(object = subset(heart, subset = condition == "BE"), group.by = "group"))

for (i in 1:length(heart.list)) {
  heart.list[[i]]@DB <- CellChatDB.human
}

future::plan("multiprocess", workers = 5)
heart.cc.list <- lapply(X = heart.list, FUN = function(x) {
  x <- subsetData(x)
  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)
  x <- computeCommunProb(x, population.size = FALSE, type =  "truncatedMean", trim = 0.1)
  x <- filterCommunication(x, min.cells = 10)
  x <- computeCommunProbPathway(x)
  x <- aggregateNet(x)
})

for (i in 1:length(heart.cc.list)) {
  heart.cc.list[[i]] <- netAnalysis_computeCentrality( heart.cc.list[[i]], slot.name = "netP")
}

cellchat <- mergeCellChat(heart.cc.list, add.names = names(heart.cc.list))
save(heart.cc.list, cellchat, file = "heart.rda")
