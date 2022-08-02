

##
##module load R
##R3.5
library(monocle)
library(Seurat)
library(ggplot2)
library(Matrix)
library(SeuratObject)
library(patchwork)

##Notes:
##The same program can be applied to other organs by inputing different seurat object.


###
outLoc = paste0("/home/ml724/project/pig")
colorRef = c("red", "deepskyblue", "seagreen", "orange", "black");

########analysis per cell type
gexpr <- NULL
gexpr = readRDS("/gpfs/ycga/project/sestan/sestanlabShare_from_Jay/OrganEx/data/OE.liver.all.seurat.filter.final.rds")

####
#unitRef = c("ALL", "BE_ECMO_H0", "BE_ECMO_H1", "BE_ECMO_H7")
unitRef = c("BE_ECMO_H0", "BE_ECMO_H1", "BE_ECMO_H7")
for(k in 1:length(unitRef)){
	eachUnit = unitRef[k];
	cat(eachUnit, "\n");
	
	###choose unit
	if(eachUnit == "ALL") condRef = c("BE", "ecmo", "h0", "h1", "h7");
	if(eachUnit == "BE_ECMO_H0") condRef = c("BE", "ecmo", "h0");
	if(eachUnit == "BE_ECMO_H1") condRef = c("BE", "ecmo", "h1");
	if(eachUnit == "BE_ECMO_H7") condRef = c("BE", "ecmo", "h7");
	
	###subset samples
	colorRef.use = alpha(colorRef[1:length(condRef)], 0.5);
	gexpr.use = subset(gexpr, subset = condition %in% condRef)
	meta = gexpr.use@meta.data

	####downsample
	Idents(gexpr.use) = "group";
	gexpr.use = subset(gexpr.use, downsample=5000);
	meta = gexpr.use@meta.data;
	
	###
	figFile = paste0(outLoc,"/pseudotime/", "pigLiver.", eachUnit, ".pseudotime.monocle2.pdf");
	pdf(figFile,10,10)
	par(omi = c(0.1, 0.1, 0.1, 0.1))
	
	###
	selectType = c("Hepatocytes")
		
	for(i in 1:length(selectType)){
		eachType = selectType[i]
		cat(eachType, "\n")
		
		###choose cell type
		gexpr2 <- NULL
		gexpr2 = subset(gexpr.use, group == eachType)
		
		###---hvg
		#gexpr2 = FindVariableFeatures(object = gexpr2, selection.method = "vst");
		#hvg = VariableFeatures(gexpr2);
		
		######-------find markers
		#topMarkers <- NULL
		#recMarkers <- NULL
		#Idents(gexpr2) = "condition"
		#
		#####
		#recMarkers <- FindAllMarkers(gexpr2, min.pct = 0.1, only.pos = F, logfc.threshold = log(1.25))
		#####----
		#idx.sig = which(as.numeric(recMarkers[,5]) < 0.01);
		#topMarkers = rownames(recMarkers)[idx.sig];
		#topMarkers = unique(recMarkers$gene)
		#cat(length(topMarkers), "\n");
				
		
		#####-----find markers
		topGenes <- NULL;
		for(j in 1:length(condRef)){
			recMarkers <- NULL;
			recMarkers <- FindMarkers(gexpr2, only.pos= F, ident.1 = condRef[j], group.by = "condition",
								min.pct = 0.1, logfc.threshold = log(1.25), max.cells.per.ident = 1000);									
			####----
			idx.sig = which(as.numeric(recMarkers[,5]) < 0.01);
			if(length(idx.sig) > 250) idx.sig = idx.sig[1:250];
			eachTop = rownames(recMarkers)[idx.sig];
			topGenes = c(topGenes, eachTop);
		}
		###
		countTop = table(topGenes);
		topMarkers = names(countTop)[which(countTop < 2)];
		#topMarkers = names(countTop);
		cat("Number of markers used: ", length(topMarkers), "\n");
		
		######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%run monocle
		##----get matrix
		meta2 = gexpr2@meta.data
		gexprMat2 = GetAssayData(gexpr2, slot="counts")
		
		
		####pheno type
		pd = new("AnnotatedDataFrame", data = meta2)
		
		###feature data
		geneSymbol = rownames(gexprMat2)
		geneAnnot = data.frame(gene_short_name=geneSymbol)
		rownames(geneAnnot) = geneSymbol
		fd = new("AnnotatedDataFrame", data = geneAnnot)
		
		#####-----expression matrix
		cds = newCellDataSet(gexprMat2, phenoData = pd, featureData = fd,
							lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
		
		####----Estimate size factors and dispersions
		cds <- estimateSizeFactors(cds)
		cds <- estimateDispersions(cds)
			

		######-----another option, using find markers foud dex genes
		cds_subset <- setOrderingFilter(cds, ordering_genes = topMarkers)
		p0 <- plot_ordering_genes(cds_subset)
		p0 <- p0 + labs(title=paste0(eachType, "_use_markers"))
		#print(p0)
		
		cds_subset <- reduceDimension(cds_subset, max_components = 2, method = 'DDRTree',norm_method = 'log')
		
		cds_subset <- orderCells(cds_subset)
		p1 <- plot_cell_trajectory(cds_subset, color_by = "condition", theta=2, cell_size=1)
		p1 <- p1 + scale_color_manual(breaks = condRef, values=colorRef.use)
		p1 <- p1 + labs(title= eachType)	
		p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=10)))
		
		
		#####-----pseudotime
		p2 <- plot_cell_trajectory(cds_subset, color_by = "Pseudotime")
		p2 <- p2 + labs(title = eachType)
		
		
		#####-----orig.ident
		p3 <- plot_cell_trajectory(cds_subset, color_by = "orig.ident", theta=2, cell_size=1)
		p3 <- p3 + labs(title = eachType)
		p3 <- p3 + guides(colour = guide_legend(override.aes = list(size=10)))
		
		#####-----State
		p4 <- plot_cell_trajectory(cds_subset, color_by = "State", theta=2, cell_size=1)
		p4 <- p4 + labs(title = eachType)
		p4 <- p4 + guides(colour = guide_legend(override.aes = list(size=10)))
		
		####----plots
		plts <- wrap_plots(list(p1, p2, p3,p4), ncol=2)
		print(plts)
	}
	dev.off()
}

