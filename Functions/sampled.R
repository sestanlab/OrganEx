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

# hipo
hipo <- readRDS(paste0(inputdir,'OE.HIP.all.seurat.filter.final.rds'))
hipo@meta.data$cdt_gp <- paste0(as.character(hipo@meta.data$condition), "_", as.character(hipo@meta.data$group))
hipo.list <- SplitObject(hipo, split.by = "cdt_gp")

set.seed(0)
hipo.list.sampled <- lapply(X = hipo.list, N = 1000, FUN = function(x, N) {
  if(ncol(x) > N){
    x <- x[,sample(ncol(x), N, replace = FALSE)] 
  } else {
    x
  }
  return(x)
})

hipo.sampled <- Reduce(merge,hipo.list.sampled)
saveRDS(hipo.sampled, "hipo.sampled.rds")

# kidney
kidney <- readRDS(paste0(inputdir,'OE.kidney.all.seurat.filter.final.rds'))
kidney@meta.data$cdt_gp <- paste0(as.character(kidney@meta.data$condition), "_", as.character(kidney@meta.data$group))
kidney.list <- SplitObject(kidney, split.by = "cdt_gp")

set.seed(0)
kidney.list.sampled <- lapply(X = kidney.list, N = 1000, FUN = function(x, N) {
  if(ncol(x) > N){
    x <- x[,sample(ncol(x), N, replace = FALSE)] 
  } else {
    x
  }
  return(x)
})

kidney.sampled <- Reduce(merge,kidney.list.sampled)
saveRDS(kidney.sampled, "kidney.sampled.rds")

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