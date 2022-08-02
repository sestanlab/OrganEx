# liver
liver <- readRDS(paste0(inputdir,'OE.liver.all.seurat.filter.final.rds'))
liver@meta.data$cdt_gp <- paste0(as.character(liver@meta.data$condition), "_", as.character(liver@meta.data$group))
liver.list <- SplitObject(liver, split.by = "cdt_gp")

set.seed(0)
liver.list.sampled <- lapply(X = liver.list, N = 1000, FUN = function(x, N) {
  if(ncol(x) > N){
    x <- x[,sample(ncol(x), N, replace = FALSE)] 
  } else {
    x
  }
  return(x)
})

liver.sampled <- Reduce(merge,liver.list.sampled)
saveRDS(liver.sampled, "liver.sampled.rds")

liver <- readRDS(paste0(inputdir, 'liver.sampled.rds'))
liver.list <- list(h0= createCellChat(object = subset(liver, subset = condition == "h0"), group.by = "group"),
                  h1= createCellChat(object = subset(liver, subset = condition == "h1"), group.by = "group"),
                  h7= createCellChat(object = subset(liver, subset = condition == "h7"), group.by = "group"),
                  ecmo= createCellChat(object = subset(liver, subset = condition == "ecmo"), group.by = "group"),
                  be=createCellChat(object = subset(liver, subset = condition == "BE"), group.by = "group"))

for (i in 1:length(liver.list)) {
  liver.list[[i]]@DB <- CellChatDB.human
}

future::plan("multiprocess", workers = 5)
liver.cc.list <- lapply(X = liver.list, FUN = function(x) {
  x <- subsetData(x)
  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)
  x <- computeCommunProb(x, population.size = FALSE, type =  "truncatedMean", trim = 0.1)
  x <- filterCommunication(x, min.cells = 10)
  x <- computeCommunProbPathway(x)
  x <- aggregateNet(x)
})

for (i in 1:length(liver.cc.list)) {
  liver.cc.list[[i]] <- netAnalysis_computeCentrality( liver.cc.list[[i]], slot.name = "netP")
}

cellchat <- mergeCellChat(liver.cc.list, add.names = names(liver.cc.list))
save(liver.cc.list, cellchat, file = "liver.rda")
