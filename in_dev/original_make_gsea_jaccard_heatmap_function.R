# Initiate -------
.libPaths('/Users/zhangy68/R/R_4.4.2')
packages_to_load <- c("dplyr", "tidyr", "fgsea", "tibble", "grid", "circlize", "dendextend")
invisible(suppressPackageStartupMessages({lapply(packages_to_load, library, character.only = TRUE)}))
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
               "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
               "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))
# Functions -------
split_genes <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character(0))
  genes <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]
  unique(genes)
}

jaccard_dist <- function(A, B) {
  A <- unique(A); B <- unique(B)
  if (length(A) == 0 && length(B) == 0) return(0)
  u <- union(A, B)
  if (length(u) == 0) return(0)
  1 - (length(intersect(A, B)) / length(u))
}


# Set input ------
# set the thresholds to only show pathways that are meaningful
threshold_positive = 2
threshold_negative = -2

# set input pathway table
lst = readRDS('RDS_files/activity.session_10.DEG_and_GSEA.scLevel.highVsLow.monocyteSubLevel1.wilcox.minpct10%.fullList.rds')
lst = lapply(seq_along(lst),function(idx){
  dfname = names(lst)[idx]
  df <- lst[[idx]] %>% mutate(Contrast = dfname)
})
# writing it in the cluster profiler output style
df.full <- lst %>%
  bind_rows %>%
  transmute(
    Contrast = Contrast,
    Database = geneSetType,
    Pathway = pathway,
    ID = pathway,
    Label = pathway,
    setSize = size,
    enrichmentScore = ES,
    NES = NES,
    pvalue = pval,
    p.adjust = padj,
    qvalue = padj,
    rank = NA,
    leading_edge = NA,
    core_enrichment = vapply(leadingEdge, function(x) paste(x, collapse = ","),character(1))
  )

df <- df.full %>% filter(p.adjust < 0.05) %>% filter(NES > threshold_positive | NES < threshold_negative)



# Read gene set gmt files --------
path.resources = '/Users/zhangy68/OneDrive - National Institutes of Health/Resources/msigdb_v2025.1.Hs_GMTs'
geneSet.reactome <- gmtPathways(file.path(path.resources, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"))
geneSet.go <- gmtPathways(file.path(path.resources, "c5.go.bp.v2025.1.Hs.symbols.gmt"))
geneSet.hallmark <- gmtPathways(file.path(path.resources, "h.all.v2025.1.Hs.symbols.gmt"))
geneSet.kegg <- gmtPathways(file.path(path.resources, "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"))
# inside geneSets, each gene set is a named vector of genes
geneSets <- c(geneSet.reactome,geneSet.go,geneSet.hallmark,geneSet.kegg)

# Attach the full gene lists of each pathway to my df as df$full_gene_list -----

# Sanity check ------
## Assumptions:
## 1) For each (Pathway, Contrast), there is at most ONE row (no duplicates).
## 2) Each Pathway has ONE identical gene-set in full_gene_list across contrasts.
##    Use full_gene_list (comma-separated genes) for Jaccard distance.

# Safety check: which pathways are missing from GMT?
missing_pathways <- setdiff(unique(df$Pathway), names(geneSets))
if (length(missing_pathways) > 0) {warning("The following pathways were not found in GMT files:\n",paste(missing_pathways, collapse = "\n"))}

df <- df %>%
  mutate(
    full_gene_list = geneSets[Pathway]
  )

# Each pathway has exactly one gene set
stopifnot(
  df %>%
    distinct(Pathway, full_gene_list) %>%
    count(Pathway) %>%
    pull(n) %>%
    all(. == 1)
)

# Each (Pathway, Contrast) appears at most once
stopifnot(
  df %>%
    count(Pathway, Contrast) %>%
    pull(n) %>%
    all(. == 1)
)



# Jaccard row clustering (by full_gene_list) + NES heatmap ------
# df must have: Contrast, Pathway, NES, full_gene_list
stopifnot(all(c("Contrast", "Pathway", "NES", "full_gene_list") %in% colnames(df)))


## 1) Pathway-level gene sets (from full_gene_list) ------
pathway_sets <- df %>%
  distinct(Pathway, full_gene_list) %>%   # safe because of Assumption 2
  mutate(genes = full_gene_list)

gene_sets <- pathway_sets$genes
names(gene_sets) <- pathway_sets$Pathway
pathways <- names(gene_sets)


## 2) Jaccard distance among pathways + hierarchical clustering ------
n <- length(pathways)
dmat <- matrix(0, nrow = n, ncol = n, dimnames = list(pathways, pathways))

# total <- n * (n + 1) / 2   # It counts how many ordered index pairs (i,j) you can form from 1..n when j≥i
# done <- 0L
# step <- max(1L, floor(total / 100))  # update ~1% increments

for (i in seq_len(n)) {
  for (j in i:n) {
    d <- jaccard_dist(gene_sets[[i]], gene_sets[[j]])
    dmat[i, j] <- d
    dmat[j, i] <- d

    # # print out the progress on distance matrix calculation
    # done <- done + 1L
    # if (done %% step == 0L || done == total) {
    #   cat(sprintf("\r%d / %d (%.1f%%) completed",
    #               done, total, 100 * done / total))
    #   flush.console()
    # }
  }
}
cat("\n")
hc <- hclust(as.dist(dmat), method = "average")

# color branches of the tree under a distance threshold
hcut = NA # FLAG for not cutting the tree
hcut <- 0.8  # distance threshold on the dendrogram height scale (Jaccard distance)
dend <- as.dendrogram(hc)
# Cut the dendrogram at height = hcut, then color branches by resulting clusters
k <- length(unique(cutree(hc, h = hcut)))
dend_col <- color_branches(dend, k = k,col = col21[1:k])
# (Optional) make branch lines thicker for visibility
dend_col <- set(dend_col, "branches_lwd", 1.2)

# # Optional: dendrogram alone
pdf("full_gene_list_dendrogram.pdf", width = 14, height = 15)  # adjust size as needed
plot(hc, main = "Pathway clustering (Jaccard distance of full_gene_list)", xlab = "", sub = "", cex = 0.7)
dev.off()
## 3) NES matrix (rows=Pathway, cols=Contrast), missing => 0 ------
nes_mat <- df %>%
  select(Pathway, Contrast, NES) %>%
  pivot_wider(names_from = Contrast, values_from = NES, values_fill = 0) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

# Order rows by Jaccard dendrogram (important step! if missing, the dendrogram label does not match to heatmap)
nes_mat = nes_mat[rownames(dmat),]
stopifnot(all(rownames(dmat) == rownames(nes_mat)))
stopifnot(all(rownames(dmat) == hc$labels))

# Ensure heatmap includes all clustered pathways (add all-zero rows if needed)
stopifnot(setequal(rownames(nes_mat),pathways))
# missing_in_nes <- setdiff(pathways, rownames(nes_mat))
# if (length(missing_in_nes) > 0) {
#   add <- matrix(0, nrow = length(missing_in_nes), ncol = ncol(nes_mat),
#                 dimnames = list(missing_in_nes, colnames(nes_mat)))
#   nes_mat <- rbind(nes_mat, add)
# }

#
## 4) Heatmap with your row clustering ------
#
# ht <- Heatmap(
#   nes_mat,
#   name = "NES",
#
#   # 1) Smaller pathway labels
#   row_names_gp = gpar(fontsize = 8),   # try 5–7
#   row_names_max_width = unit(20, "cm"),# allow long names to extend (avoid too much squeezing)
#
#   # 2) Make dendrogram big and visible (about half the plot)
#   cluster_rows = as.dendrogram(hc),
#   show_row_dend = TRUE,
#   row_dend_width = unit(12, "cm"),     # increase/decrease (e.g., 8–12 cm)
#
#   # 3) Heatmap has only 3 columns -> make it narrow
#   cluster_columns = T,
#   width = unit(2.5, "cm"),             # adjust (e.g., 2–4 cm)
#
#   row_names_side = "left",
#   column_names_side = "top",
#
#   # Optional: make column names readable
#   column_names_gp = gpar(fontsize = 10)
# )
#
# draw(ht)


# Discrete NES mapping

nes_disc <- nes_mat
nes_disc[nes_disc >  threshold_positive] <-  1
nes_disc[nes_disc <  threshold_negative] <- -1

col_fun_disc <- colorRamp2(
  c(-1, 0, 1),
  c("blue", "white", "red")
)

stopifnot(all(rownames(dmat) == rownames(nes_disc)))

# assign the proper dendrogram based on the hcut value
row_clust <- if (is.na(hcut)) as.dendrogram(hc) else dend_col


ht <- Heatmap(
  nes_disc,
  name = "NES",

  col = col_fun_disc,

  rect_gp = gpar(col = "black", lwd = 0.4),
  # 1) Smaller pathway labels
  row_names_gp = gpar(fontsize = 6,fontface = 'bold'),
  row_names_max_width = unit(20, "cm"),

  # 2) Make dendrogram big and visible
  cluster_rows = row_clust,
  row_dend_reorder = FALSE,
  show_row_dend = TRUE,
  row_dend_width = unit(12, "cm"),
  row_split = ifelse(is.na(hcut),NULL,k),

  # 3) Heatmap has only 3 columns -> make it narrow
  cluster_columns = F,
  width = unit(2.5, "cm"),

  row_names_side = "left",
  column_names_side = "top",

  column_names_gp = gpar(fontsize = 10),

  # Make legend explicit for discrete bins
  heatmap_legend_param = list(
    at = c(-1, 0, 1),
    labels = c(
      paste0("< ", threshold_negative),
      paste0("[", threshold_negative, ", ", threshold_positive, "]"),
      paste0("> ", threshold_positive)
    )
  )
)

pdf("NES_heatmap.pdf", width = 14, height = 15)  # adjust size as needed

draw(ht)
dev.off()
