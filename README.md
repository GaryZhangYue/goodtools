
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!--  use devtools::build_readme() to update -->

# goodtools

<!-- badges: start -->
<!-- badges: end -->

The goal of goodtools is to …

## Installation

You can install the development version of goodtools like so:

``` r
remotes::install_github("GaryZhangYue/goodtools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(goodtools)
# res = make_gsea_jaccard_heatmap(df.full = df.full,
#                           padj_cutoff = padj_cutoff,
#                           threshold_positive = threshold_positive, threshold_negative = threshold_negative,
#                           hclust_method = 'average',
#                           hcut = hcut,
#                           branch_palette = col21,
#                           branches_lwd = 2,
#                           row_split = row_split,
#                           cluster_columns = F,
#                           heatmap_width_cm = 2.5,   
#                           row_dend_width_cm = 6,
#                           row_label_fontsize = 6,
#                           row_label_bold = TRUE,
#                           row_names_max_width_cm = 20,
#                           cell_border_col = "black",
#                           cell_border_lwd = 0.4,
#                           pdf_heatmap = paste0("NES_heatmap.monocyteL1.padj", padj_cutoff, ".pos", threshold_positive, 
#                                                ".neg", threshold_negative, ".rowsplit_", row_split, ".hcut", hcut, ".goodtools.pdf"),
#                           pdf_dendrogram = paste0("dendrogram.monocyteL1.padj", padj_cutoff, ".pos", threshold_positive, 
#                                                   ".neg", threshold_negative, ".rowsplit_", row_split, ".hcut", hcut, ".goodtools.pdf"),
#                           pdf_width = 12,
#                           pdf_height = 9,
#                           plot_dendrogram = TRUE,
#                           show_progress = FALSE)
```
