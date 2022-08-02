library(RColorBrewer);
library(circlize);
library(ComplexHeatmap);
library(GetoptLong);
library(ggpubr);
library(WGCNA);
library(Seurat);
library(GenomicRanges);
library(GenomicFeatures);
library(GenomicAlignments);
library(biomaRt);
library(AnnotationDbi);
library(org.Ss.eg.db);
library(topGO);
library(pheatmap);
library(igraph);
library(tidyverse);
options(stringsAsFactors=FALSE);
source("genecluster.function.R");
## use args to have tissue and cell types
tissue="heart";
celltype="Cardiomyocytes";
## dir for seurat files
sratfiledir="/gpfs/gibbs/pi/sestan.ycga/sestanlabShare_from_Jay/OrganEx/data/Final_120321/";
# show files
#dir("/gpfs/gibbs/pi/sestan.ycga/sestanlabShare_from_Jay/OrganEx/data/Final_120321/")
#[1] "Final_12032021_OE.heart.highReso.rds"
#[2] "Final_12032021_OE.HIP.highReso.rds"
#[3] "Final_12032021_OE.kidney.highReso.rds"
#[4] "Final_12032021_OE.liver.highReso.rds"
## start with HIP neurons
tissue_srat=readRDS(paste0(sratfiledir,"Final_12032021_OE.",tissue,".highReso.rds"));
tissue_celltype=subset(tissue_srat,subset=group==celltype);
Idents(tissue_celltype)="condition";
conditions = unique(tissue_celltype@meta.data$condition);
BEvsAll=data.frame();
for (i in 1:(length(conditions)-1)){
    for (j in (i+1):length(conditions)){
        thisBEvs=FindMarkers(tissue_celltype, ident.1 = conditions[i], ident.2 = conditions[j], max.cells.per.ident = 5000);
#        thisBEvs$cutFC1=(thisBEvs$pct.1+0.01)/(thisBEvs$pct.2+0.01);
#        thisBEvs$cutFC2=(thisBEvs$pct.2+0.01)/(thisBEvs$pct.1+0.01);
        ## for heart
        thisBEvs=thisBEvs[which(thisBEvs$p_val_adj < 0.05),];
        thisBEvs=thisBEvs[abs(thisBEvs$avg_log2FC) > 0.75,];
        thisBEvs$gene=rownames(thisBEvs);
        BEvsAll=rbind(BEvsAll,thisBEvs);
    }
}
save(BEvsAll,file=paste0(tissue,"_",celltype,"_selectedGenes.Rdata"));
## unique genes
seletedExpr=AverageExpression(tissue_celltype,features =  unique(BEvsAll$gene), group.by = "samplename");
data=seletedExpr$RNA;
minClusterSize = 50
sensitivity = 2
tree.method = "ward.D2"
cor_method = "p"
para <- list(minClusterSize = minClusterSize, tree.method = tree.method, sensitivity = sensitivity, cor_method = cor_method)
data_use <- t(data)
dissTOM <- 1-cor(data_use, method = cor_method)
geneTree <- dissTOM %>% as.dist() %>% hclust(method = tree.method)

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = sensitivity, pamRespectsDendro = FALSE,minClusterSize = minClusterSize);
new_mods <- dynamicMods
names(new_mods)[names(new_mods) == ""] <- "0"
mes <- moduleEigengenes(data_use, colors = paste0("_", new_mods))$eigengenes
mes <- mes %>% as.matrix() %>% scale() %>% MinMax(., min = -2.5, max = 2.5) %>% as.data.frame(., check.names = FALSE)
#mes$condition = sapply(rownames(mes),function(x) unlist(strsplit(x,"_"))[2]);
## Get the gene annotation
gene_label <- setNames(paste0("ME_", dynamicMods), colnames(data_use))
gene_label[gene_label == "ME_"] <- "ME_0"
mm <- as.data.frame(cor(data_use, mes, use = "p"));
gene_meta <- Anno_gene(gene_label = gene_label, mm = mm)
mes$condition = sapply(rownames(mes),function(x) unlist(strsplit(x,"_"))[2]);
modules = unique(gene_meta$module);
write.csv(mes,"heart_Cardiomyocytes_modulesgenes.csv");
## plot network
pdf(paste0(tissue,"_",celltype,"_network.pdf"));
weights=matrix(1,nrow=dim(dissTOM)[1],ncol=dim(dissTOM)[2]);
for (i in 1:length(modules)){
    weights[which(gene_meta$module==modules[i]),which(gene_meta$module!=modules[i])]=0;
}
netmat=matrix(0,nrow=dim(dissTOM)[1],ncol=dim(dissTOM)[2]);
colnames(netmat)=colnames(dissTOM);
rownames(netmat)=rownames(dissTOM);

netmat[which(dissTOM < 0.50 & weights==1, arr.ind=T)]=1;
netmat[which(dissTOM < 0.20 & weights==0, arr.ind=T)]=1;

diag(netmat)=0;
modulecol = brewer.pal(length(modules), "Set3");
gene_meta$col=modulecol[as.factor(gene_meta$module)];

net <- graph_from_incidence_matrix(netmat, directed=F);
V(net)$color=gene_meta$col;
V(net)$size=3;
set.seed(20211207);
lo=layout_with_fr(net);
#Isolated = which(degree(net)==0);
#net = delete.vertices(net, Isolated);
#lo = lo[-Isolated,];
plot(net,vertex.label=NA,layout=lo);
legend(x = "bottomright", legend=modules, cex=1.25 * 1, col=modulecol[as.factor(modules)], pch=rep(20, 4), title = "Modules")
dev.off();

## plot
pdf(file=paste0(tissue,"_",celltype,"_Modules_TopGO_genes.pdf"));
## find protein coding promoters and non-protein coding promoters
#mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl");
#hgncdata=getBM(attributes =c("chromosome_name","transcription_start_site","strand","entrezgene_id","gene_biotype","hgnc_symbol"),mart=mart);
#save(hgncdata,file="hgncdata_Sscrofa10.2.Rdata");
load("../hgncdata_Sscrofa10.2.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
GOmoduledata=data.frame();
for (i in 1:length(modules)){
    df=mes[,c(which(colnames(mes)==modules[i]),dim(mes)[2])];
    colnames(df)[1]="x";
    df.summary <- df %>% group_by(condition) %>% summarize(ymin = mean(x)-sd(x)/sqrt(length(x)),ymax = mean(x)+sd(x)/sqrt(length(x)),ymean = mean(x));
    df.summary$condition=factor(df.summary$condition,levels=c("h0","h1","h7","ecmo","BE"));
    g=ggplot(df.summary, aes(x = condition, y = ymean)) + geom_point(size = 2) + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + ggtitle(modules[i]) + ylim(-3,3) + geom_line(group = 1);
    print(g);
    plotmat=seletedExpr$RNA[which(gene_meta$module==modules[i]),];

    targetgenes=gene_meta$gene[which(gene_meta$module == modules[i])];
    ##All genes
    geneNames= unique(codgene$entrezgene_id);
    ##top100 upregulated gene
    output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% targetgenes)]);
    myInterestingGenes = output;
    geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
    names(geneList) <- geneNames;
    GOdata <- new("topGOdata",ontology = "BP",description = "Module genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Ss.eg.db");
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher");
    UpallRes <- GenTable(GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
    UpallRes$module=modules[i];
    GOmoduledata=rbind(GOmoduledata,UpallRes);
    ## plot
    barplot(-log10(as.numeric(UpallRes$classicFisher)), main=paste0(tissue,"_",modules[i]," genes"), horiz=TRUE, names.arg=UpallRes$Term, las=2);
    abline(v=-log10(0.01),col="red");
}
write.csv(GOmoduledata,file="heart_Cardiomyocytes_GO_modules.csv");
dev.off();
