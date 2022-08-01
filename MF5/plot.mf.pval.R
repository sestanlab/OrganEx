library(dplyr)
library(tibble)
library(ggplot2)


rawdata <- read.table("./load_files/CelltypePval.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
			tidyr::gather(., "organ", "p_val", c("Heart", "Hippocampus", "Kidney", "Liver")) %>%
			mutate(p_val = -log10(p_val))
rawdata$p_val[rawdata$p_val > 4] <- 4



orgs <- c("Heart", "Hippocampus", "Kidney", "Liver")
plist <- lapply(orgs, function(xx) {
	subdata <- rawdata %>%
				filter(organ == xx) %>%
				mutate(condition = factor(condition, levels = rev(c("h0", "h1", "h7", "ecmo", "BE")))) %>%
				mutate(color = ifelse(p_val >= -log10(0.05), "red", "lightgrey"))


	com_paths <- c("pos reg of DNA rep", "pos reg of cyto", "neg reg of apoptosis", "atp metabolic process")
	p_ord <- switch(xx,
				Heart = c(com_paths, "FA B ox", "glycol proc", "cardiac m cell AP"), 
				Hippocampus= c(com_paths, "astro panreactive", "microglia proinfl"), 
				Kidney = c(com_paths, "pct injury", "transporters"), 
				Liver = c(com_paths, "acute phase reac", "cyp"))

	subdata <- subdata %>%
				filter(pathway %in% p_ord) %>%
				mutate(orgidx = as.numeric(factor(pathway, levels = p_ord)))
	print(sum(is.na(subdata$orgidx)))


	p <- ggplot(subdata, aes_string(x = "orgidx", y = "condition", size = "p_val", color = "p_val")) +
				geom_point(shape = 16)+
				scale_radius(range = c(0, 6)) +
				scale_color_gradient2(low = "lightgrey", mid = "lightgrey", high = "red", midpoint = -log10(0.05)) + 
				##scale_color_identity() + 
				scale_x_discrete(limits = 1:7, breaks = 1:max(subdata$orgidx), labels = c(p_ord))+
				coord_fixed() + 
				theme_bw() +
				theme(panel.grid.major = element_line(size= 0.2), panel.grid.minor = element_line(size = 0.2), axis.title = element_blank())
	p
	})

pdf(paste0("./report/", "MF6_Pval_dotplot_v2.pdf"), width = 8, height = 9, useDingbats = FALSE)
patchwork::wrap_plots(plist, nrow = length(orgs), ncol = 1, guides = "collect") & theme(legend.position = "right")
dev.off()













