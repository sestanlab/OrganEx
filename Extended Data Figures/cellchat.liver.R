
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
