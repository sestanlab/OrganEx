
hipo <- readRDS(paste0(inputdir, 'hipo.sampled.rds'))
hipo.list <- list(h0= createCellChat(object = subset(hipo, subset = condition == "h0"), group.by = "group"),
                  h1= createCellChat(object = subset(hipo, subset = condition == "h1"), group.by = "group"),
                  h7= createCellChat(object = subset(hipo, subset = condition == "h7"), group.by = "group"),
                  ecmo= createCellChat(object = subset(hipo, subset = condition == "ecmo"), group.by = "group"),
                  be=createCellChat(object = subset(hipo, subset = condition == "BE"), group.by = "group"))

for (i in 1:length(hipo.list)) {
  hipo.list[[i]]@DB <- CellChatDB.human
}

future::plan("multiprocess", workers = 5)
hipo.cc.list <- lapply(X = hipo.list, FUN = function(x) {
  x <- subsetData(x)
  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)
  x <- computeCommunProb(x, population.size = FALSE, type =  "truncatedMean", trim = 0.1)
  x <- filterCommunication(x, min.cells = 10)
  x <- computeCommunProbPathway(x)
  x <- aggregateNet(x)
})

for (i in 1:length(hipo.cc.list)) {
  hipo.cc.list[[i]] <- netAnalysis_computeCentrality( hipo.cc.list[[i]], slot.name = "netP")
}

cellchat <- mergeCellChat(hipo.cc.list, add.names = names(hipo.cc.list))
save(hipo.cc.list, cellchat, file = "hipo.rda")
