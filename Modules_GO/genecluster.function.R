##library(crestree);
library(WGCNA)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(ggpubr)
plot_sp_heatmap <- function(meta, data, scale.exp = FALSE, split.by = "species", group.by = "mres", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), group.order, min_max = c(-2.5, 2.5), file_name, output_dir = outputdir, pdf_width = 8, show_rownames = TRUE, height_unit = 0.1, label_genes = NULL, font_scale = 1, row_meta = NULL, return.plot = FALSE) {
    ## Set the MinMax Values
    if (is.null(group.order)){
        group.order <- levels(as.factor(as.character(meta[, group.by])))
    }
    if (is.null(split.order)){
        split.order <- levels(as.factor(as.character(meta[, split.by])))
    }
    if (scale.exp){
        data <- data %>% t() %>% scale() %>% t()
    }
    data <- data %>% MinMax(., min = min_max[1], max = min_max[2])
    nsplit <- length(split.order)


    ## Do the plots
    plot_meta <- meta[, c(group.by, split.by)]
    sp.cols <- c("#FF420E","#4CB5F5","#89DA59","#89DA59","#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Macaque", "Marmoset"))
    cls.cols <- gg_color_hue(length(group.order)) %>% setNames(., group.order)
    anno_cols <- list(sp.cols[plot_meta[, split.by] %>% as.factor() %>% levels()], cls.cols[plot_meta[, group.by] %>% as.factor() %>% levels()]) %>%
                        setNames(., c(split.by, group.by))

    ## Reorganize the plot data
    plot_data <- lapply(split.order, function(x) {
        submeta <- plot_meta[plot_meta[, split.by] == x, ]
        td <- data[, rownames(submeta), drop = FALSE] %>% as.matrix()
        rownames(td) <- paste0(x, "|", rownames(td))
        colnames(td) <- as.character(submeta[, group.by])
        td <- td[, group.order, drop = FALSE]
        return(td)
        }) %>%
            do.call(rbind, .)
    gene_orders <- lapply(rownames(data), function(x) grep(paste0("\\|", x, "$"), rownames(plot_data))) %>%
                unlist(., use.names = FALSE)
    plot_data <- plot_data[gene_orders, ]
    row_labels <- lapply(rownames(data), function(x) c("", "", x, "")) %>% unlist(., use.names = FALSE)

    if (is.null(label_genes)){
        ngenes <- nrow(data)
        pdf_heights <- max(ceiling(ngenes * height_unit), 20)
        row_fontsize <- ceiling(12/ sqrt(ngenes)) * 3 * font_scale %>% MinMax(., min = 1.5, max = 12)
        pdf(paste0(output_dir, file_name, "_heatmap.pdf"), width = pdf_width, height = pdf_heights)
        pheatmap::pheatmap(plot_data, cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(3))(30), border_color = NA, annotation_col = plot_meta, annotation_row = NA,annotation_colors = anno_cols, show_rownames = show_rownames, labels_row=row_labels, show_colnames = FALSE, fontsize_col = 8, fontsize_row = row_fontsize, gaps_row = seq(length(split.order), nrow(plot_data), length(split.order)))
        dev.off()
    } else {
        column_ha <- HeatmapAnnotation(df = plot_meta, col = anno_cols, annotation_height = unit(c(0.01, 0.01), "in"))
        ## Get the range of the heatmap legends
        legend_limits <- c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)
        color_breaks <- legend_limits[legend_limits >= min(data[!is.na(data)]) & legend_limits <= max(data[!is.na(data)])]
        font_size <- font_scale * 5
        if (is.list(label_genes)){
            htlist <- rowAnnotation(link = anno_mark(at = which(rownames(data) %in% label_genes[[1]]),
                    labels = rownames(data)[which(rownames(data) %in% label_genes[[1]])],
                    side = "left",
                    labels_gp = gpar(fontsize = font_size), padding = unit(1, "mm"))) +
                Heatmap(data, col = colorRamp2(color_breaks, viridis(length(color_breaks))), na_col = "white",
                    name = "scaled_expr", column_title = file_name, show_column_names = FALSE, width = unit(pdf_width - 4, "in"),
                    heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks), cluster_rows = FALSE, cluster_columns = FALSE, column_split = plot_meta[, split.by], top_annotation = column_ha, row_names_gp = gpar(fontsize = 5)) +
                rowAnnotation(link = anno_mark(at = which(rownames(data) %in% label_genes[[2]]),
                    labels = rownames(data)[which(rownames(data) %in% label_genes[[2]])],
                    side = "right",
                    labels_gp = gpar(fontsize = font_size), padding = unit(1, "mm")))
        } else {
            htlist <- Heatmap(data, col = colorRamp2(color_breaks, viridis(length(color_breaks))), na_col = "white",
                name = "scaled_expr", column_title = file_name, show_column_names = FALSE, width = unit(pdf_width - 4, "in"),
                heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks), cluster_rows = FALSE, cluster_columns = FALSE, column_split = plot_meta[, split.by], top_annotation = column_ha, row_names_gp = gpar(fontsize = 5)) +
                rowAnnotation(link = anno_mark(at = which(rownames(data) %in% label_genes),
                    labels = rownames(data)[which(rownames(data) %in% label_genes)],
                    labels_gp = gpar(fontsize = font_size), padding = unit(1, "mm")))
        }

        if (return.plot) {
            return(htlist)
        } else {
            pdf_heights <- 12
            pdf(paste0(output_dir, file_name, "_heatmap.pdf"), width = pdf_width, height = pdf_heights)
            draw(htlist)
            dev.off()
        }


    }
}

Anno_gene <- function(gene_label, mm) {
    ## Build the Gene Metadata
    gene_label <- gene_label[rownames(mm)]
    gene_meta <- data.frame(gene = names(gene_label),
                    module = gene_label,
                    modulemembership = sapply(names(gene_label), function(g) round(mm[g, gene_label[g]], digits = 4)),
                    stringsAsFactors = FALSE) %>%
                    mutate(id = paste0(module, "::", gene, "::", as.character(modulemembership)))
    colnames(gene_meta)[colnames(gene_meta) == "modulemembership"] <- "mm"


    ## Get the order of genes within each module
    all_modules <- levels(as.factor(gene_label))
    gene_meta$gorder <- NA
    all_order <- c()
    for (mod in all_modules){
        idx <- which(gene_label == mod)
        cur_order <- gene_meta$mm[idx] %>% setNames(., gene_meta$gene[idx]) %>% sort(decreasing = TRUE) %>% names()
        gene_meta$gorder[idx] <- match(cur_order, gene_meta$gene)
        all_order <- c(all_order, cur_order)
    }
    gene_meta$all_order <- match(all_order, gene_meta$gene)
    rownames(gene_meta) <- gene_meta$gene
    return(gene_meta)
}




get_modules <- function(data, meta, minClusterSize = 5, sensitivity = 4, file_name, split.by = "species", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), group.by = "mres", group.order = NULL, output_dir = outputdir){
    ## Store the parameters [in a list]
    tree.method = "ward.D2"
    cor_method = "p"
    para <- list(minClusterSize = minClusterSize, tree.method = tree.method, sensitivity = sensitivity, cor_method = cor_method)


    ##Cluster the genes [Remove the na_cols first]
    data_use <- t(data)
    dissTOM <- 1-cor(data_use, method = cor_method)
    geneTree <- dissTOM %>% as.dist() %>% hclust(method = tree.method)


    ## Get the gene modules
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = sensitivity, pamRespectsDendro = FALSE,
                            minClusterSize = minClusterSize);
    new_mods <- dynamicMods
    names(new_mods)[names(new_mods) == ""] <- "0"
    mes <- moduleEigengenes(data_use, colors = paste0("_", new_mods))$eigengenes
    mes <- mes %>%
            as.matrix() %>% scale() %>%
            MinMax(., min = -2.5, max = 2.5) %>%
            as.data.frame(., check.names = FALSE)
    ##print(rownames(mes))


    ## Get the gene annotation
    gene_label <- setNames(paste0("ME_", dynamicMods), colnames(data_use))
    gene_label[gene_label == "ME_"] <- "ME_0"
    mm <- as.data.frame(cor(data_use, mes, use = "p"));
    gene_meta <- Anno_gene(gene_label = gene_label, mm = mm)



    ## Visualize the eigengenes
    mod_order <- gsub("ME_", "", colnames(mes)) %>% as.numeric() %>%
                sort() %>% paste0("ME_", .)
    plot_sp_heatmap(meta = meta, data = as.matrix(t(mes[, mod_order])), split.by = split.by, split.order = split.order, group.by = group.by, group.order = group.order, file_name = paste0(file_name, "_eigen"), output_dir = output_dir)

    ## Plot the raw expression
    plot_sp_heatmap(meta = meta, data = data[rownames(gene_meta)[gene_meta$all_order], ], scale.exp = TRUE, split.by = split.by, split.order = split.order, group.by = group.by, group.order = group.order, file_name = paste0(file_name, "_allgenes"), output_dir = output_dir)


    ## Plot the tree of hierachical clustering
    hc <- hclust(as.dist(1-cor(mes)), method = tree.method)
    pdf(paste0(output_dir, file_name, "_module_tree.pdf"), width = 5, height = 4)
    p <- ggdendro::ggdendrogram(hc) + geom_hline(yintercept = 0.25, col = "red")
    print(p)
    dev.off()


    out_res <- list(data = data, meta = meta, mes = mes[colnames(data),], mm = mm[rownames(data),], gene_meta = gene_meta[rownames(data),], hc = hc)
    return(out_res)
}
