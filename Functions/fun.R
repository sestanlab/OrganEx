
# hyper.fun.R

hp_test <- function(all_genes, marker_list, test_genes){
  total_white <- length(test_genes)
  total_black <- length(all_genes) - total_white
  sampling_size <- sapply(marker_list, length)
  capture <- sapply(marker_list, function(x) sum(x %in% test_genes))
  pval <- sapply(1:length(marker_list), function(x) {
    pp <- sum(dhyper(capture[x]:sampling_size[x], total_white, total_black, sampling_size[x]))
    return(pp)
  }) %>%
    setNames(., names(marker_list))
  return(pval)
}