## Put the metadata together
source("~/project/HIP/scripts/hip.fun.R")

seu_list <- list(Heart = readRDS(file = "../Heart/load_files/OE.heart.all.seurat.filter.final.rds"),
				Hippocampus = readRDS(file = "../Hippocampus/load_files/OE.HIP.all.seurat.filter.final.rds"),
				Kidney = readRDS("../Kidney/load_files/OE.kidney.all.seurat.filter.final.rds"),
				Liver = readRDS("../Liver/load_files/OE.liver.all.seurat.filter.final.rds"))



meta <- lapply(names(seu_list), function(organ) {
	subm <- seu_list[[organ]]@meta.data
	subm <- subm[, c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "samplename", "cluster", "group")]
	subm$organ <- organ
	return(subm)
	}) %>%
			do.call(rbind, .)
saveRDS(meta, file = paste0(inputdir, "OE.meta.all.rds"))



###----------------------------------------------------------------------------------------
## Plot number of cells by samples & quality metrics
library(lemon)


## Further add a replicate column
meta <- readRDS(file = paste0(inputdir, "OE.meta.all.rds"))

meta <- meta %>%
			mutate(donor = extract_field(samplename, 1, "_"))

plot_data <- meta %>%
				rownames_to_column("cell") %>%
				group_by(organ, samplename) %>%
				summarize(size = n()) %>%
				ungroup() %>%
				mutate(condition = extract_field(samplename, 2, "_")) %>%
				mutate(org_condi = paste0(organ, "|", condition))


ind.cells <- plot_data %>% 
				group_by(org_condi) %>%
				summarize(size = sum(size)) %>%
				ungroup()


sample_order <- lapply(c("Heart", "Hippocampus", "Kidney", "Liver"), function(organ){
	paste0(organ, "|", c("h0", "h1", "h7", "ecmo", "BE"))
	}) %>%
		unlist()


cdt_color <- c("#6650A2", "#E1BD8E", "#DFD4E9", "#FFA500", "#CD5C5C") %>%
            setNames(., c("BE", "h0", "h1", "h7", "ecmo"))

org_ncells <- plot_data %>%
				group_by(organ) %>%
				summarize(size = sum(size)) %>%
				ungroup() %>%
				mutate(y_loc = 55000) %>%
				mutate(org_condi = paste0(organ, "|", "h7"))


p <- ggplot(plot_data, aes_string(x = "org_condi", y = "size")) +
		geom_bar(aes_string(fill = "condition"), color = "black", position = position_stack(reverse = FALSE), stat = "identity", lwd = 0.3) +
		geom_text(data = ind.cells, mapping = aes_string(x = "org_condi", y = "size", label = "size"), nudge_y = 1500, angle = 45, vjust = 0.5, hjust = 0.1, size = 2.8) +
		geom_label(data = org_ncells, aes_string(x = "org_condi", y = "y_loc", label = "size"), fill = NA) + 
		coord_capped_cart(left='both') +
		scale_x_discrete(limits = sample_order) +  	
		scale_y_continuous(limits = c(0, 60000)) +
		scale_fill_manual(values = cdt_color) +
		theme_classic() + 
		RotatedAxis() + 
		labs(y = "Sample size", x = "Individual") +
		theme(axis.line=element_line(size = 0.25), axis.ticks=element_line(size = 0.25)) +
		theme(axis.text.x = element_text(size = rel(0.9)),  axis.title.x= element_blank())

pdf(paste0(outputdir, "OE.sample.size.pdf"), width = 8, height = 4)
print(p)
dev.off()



## Violin plot comparing the quality (nFeature, nCount)
qc_data <- meta %>%
			mutate(condition = extract_field(samplename, 2, "_")) %>%
			mutate(org_condi = paste0(organ, "|", condition)) %>%
			group_by(org_condi) %>%
			mutate(donoridx = paste0("s", as.numeric(as.factor(samplename)))) %>%
			ungroup()

plist <- lapply(c("nCount_RNA","nFeature_RNA"), function(x) {
		p2 <- ggplot(qc_data, aes_string(x = "org_condi", y = x, fill = "donoridx")) +
			geom_violin(position='dodge', scale = "width", size = 0.05, adjust = 1.5, trim =TRUE, alpha = 0.8) + 
			coord_capped_cart(top='both', left='both') +
			scale_x_discrete(limits = sample_order) +  	
			theme_classic() + 
			RotatedAxis() + 
			labs(y = switch(x, nFeature_RNA = "nGenes", nCount_RNA = "nUMIs")) +
			theme(panel.border=element_blank(), axis.line=element_line(size = 0.2), axis.ticks = element_line(size = 0.2),axis.text.x = element_text(size = rel(0.9)), axis.text.y = element_text(size = rel(0.9)), axis.title.x = element_blank(), axis.title.y = element_text(size = rel(.8)), legend.position = "none")
		p2
		})
plist[[1]] <- plist[[1]] +
			theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

pdf(paste0(outputdir, "OE.Sample.quality.pdf"), width = 8, height = 7)##, units = "in",  res = 300)
patchwork::wrap_plots(plist, nrow = 2, ncol = 1)
dev.off()






