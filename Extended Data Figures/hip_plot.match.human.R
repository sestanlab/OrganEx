source("~/project/HIP/scripts/hip.fun.R")
library(ggdendro)



##--------------------------------------------------------------------------------
## heatmap visualization of dataset comparisons 
hipec <- readRDS(file = paste0("~/project/HIP/data/", "HIPEC.final.seu.01172021.rds"))
hvg_human <- rownames(hipec$pca@feature.loadings)
exp_human <- rownames(hipec)[Matrix::rowSums(hipec$RNA@counts!=0) >= 20]


avgs_human <- readRDS(file = paste0("~/project/HIP/MF1_overview/load_files/", "Avg.HIPEC.subtypes.rds"))
cls_order <- c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","SUB proximal ROBO1 SEMA3E","SUB distal FN1 NTNG1","EC L2 CUX2 PDGFD","EC L2 CUX2 LAMA3","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L6 THEMIS CDH13","EC L5 RORB TPBG","EC L6 THEMIS RGS12","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2","CR RELN NDNF","InN MEIS2 SHISAL2B","InN SST NPY","InN SST ADAMTS12","InN SST EPB41L4A","InN SST OTOF","InN PVALB PLEKHH2","InN LHX6 AC008415.1","InN PVALB MEPE","InN PVALB PLCL1","InN VIP NOX4","InN VIP SCTR","InN VIP ABI3BP","InN VIP SCML4","InN VIP CHRNA2","InN VIP PENK","InN NR2F2 PTPRK","InN NR2F2 MIR4300HG","InN NR2F2 SLC17A8","InN NR2F2 ANO2","InN NR2F2 DDR2","InN LAMP5 CHST9","InN LAMP5 KIT","InN LAMP5 NMBR","Astro AQP4 CHRDL1","Astro AQP4 GFAP","OPC PDGFRA EGR1","OPC PDGFRA GRIA4","COP GPR17 ADAM33","Oligo CPXM2 KANK4","Oligo OPALIN LAMA2","Oligo OPALIN LINC01098","Oligo OPALIN SLC5A11","Micro C1QB CD83","Micro C1QB P2RY12","Macro F13A1 COLEC12","Myeloid LSP1 LYZ","T SKAP1 CD247","aEndo DKK2 FBLN5","Endo CLDN5 VWF","PC CLDN5 ABCC9","vSMC ABCC9 P2RY14","aSMC ACTA2 TAGLN","VLMC COL1A1 COL1A2")    


pig <- readRDS(file = paste0(inputdir, "OE.HIP.all.seurat.filter.final.highReso.rds"))
pig@meta.data$cluster <- gsub("GC_neuroblast|GC_mature", "GC", pig@meta.data$hres) %>%
                        gsub("CA2-3", "CA2-4", .) %>%
                        gsub("EC_up|EC_deep", "EC", .) %>%
                        gsub("InN_MGE|InN_CGE", "InN", .)
hvg_pig <- subset(pig, condition == 'h0') %>%
            SplitObject(., split.by = "samplename") %>%
            lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
            SelectIntegrationFeatures(., nfeatures = 2000)


Idents(pig) <- "cluster"
avg_pig <- log(AverageExpression(pig)$RNA + 1)
exp_pig <- rownames(pig)[Matrix::rowSums(pig$RNA@counts!=0) >= 15]


pig_order <- c("GC","Mossy_cells","CA2-4", "CA1_Sub","EC","InN","Astro","OPC","Oligo","Micro","Vas")


gene_use <- intersect(hvg_pig, exp_human)
cor_mat <- cor(avg_pig[gene_use, ], avgs_human[gene_use, ], method = "p")


pdf(paste0(outputdir, "SF1.HIP.human-pig.heatmap v2.pdf"), width = 7, height = 5)
pheatmap::pheatmap(cor_mat[pig_order, cls_order], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 5, fontsize_row = 8)
dev.off()






##--------------------------------------------------------------------------------
## Plot UMAP
pig <- readRDS(file = paste0(inputdir, "OE.HIP.all.seurat.filter.final.highReso.rds"))
gp_colors <- c("#4472C4", "#70AD47", "#E0BC8D", "#FFC000", "#C83200") %>% 
          setNames(., c(c("Neuron", "Astro", "OPC Oligo", "Micro", "Vas")))
DimFig(pig, group.by = "group", cols = gp_colors, file_name = "SF1.HIP.lowreso.UMAP", plot.scale = 0.9)


cls_cols <- setNames(c("#31028e", "#5b178c", "#8d48af", "#633ea3", "#8100dd", "#7a17c1", "#E1BD8E", "#deccff", "#d9aff7", "#FFA500", "#CD5C5C"), 
	c("GC","Mossy_cells","CA2-4", "CA1_Sub","EC","InN","Astro","OPC","Oligo","Micro","Vas"))
DimFig(pig, group.by = "cluster", cols = cls_cols, file_name = "SF1.HIP.highreso.UMAP", plot.scale = 0.9)















