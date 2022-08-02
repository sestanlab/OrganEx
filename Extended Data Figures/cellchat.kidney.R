
kidney <- readRDS(paste0(inputdir, 'kidney.sampled.rds'))
kidney.list <- list(h0= createCellChat(object = subset(kidney, subset = condition == "h0"), group.by = "group"),
                  h1= createCellChat(object = subset(kidney, subset = condition == "h1"), group.by = "group"),
                  h7= createCellChat(object = subset(kidney, subset = condition == "h7"), group.by = "group"),
                  ecmo= createCellChat(object = subset(kidney, subset = condition == "ecmo"), group.by = "group"),
                  be=createCellChat(object = subset(kidney, subset = condition == "BE"), group.by = "group"))

for (i in 1:length(kidney.list)) {
  kidney.list[[i]]@DB <- CellChatDB.human
}

future::plan("multiprocess", workers = 5)
kidney.cc.list <- lapply(X = kidney.list, FUN = function(x) {
  x <- subsetData(x)
  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)
  x <- computeCommunProb(x, population.size = FALSE, type =  "truncatedMean", trim = 0.1)
  x <- filterCommunication(x, min.cells = 10)
  x <- computeCommunProbPathway(x)
  x <- aggregateNet(x)
})

for (i in 1:length(kidney.cc.list)) {
  kidney.cc.list[[i]] <- netAnalysis_computeCentrality( kidney.cc.list[[i]], slot.name = "netP")
}

cellchat <- mergeCellChat(kidney.cc.list, add.names = names(kidney.cc.list))
save(kidney.cc.list, cellchat, file = "kidney.rda")
