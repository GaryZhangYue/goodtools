#' Jaccard-based pathway clustering and discrete NES heatmap
#'
#' Perform hierarchical clustering of pathways using Jaccard distance on
#' full gene sets, and visualize discrete normalized enrichment scores (NES)
#' across contrasts as a heatmap with an aligned dendrogram.
#'
#' This function assumes that `df.full` has already been prepared and
#' contains a `full_gene_list` column (one identical gene set per pathway
#' across contrasts). No GMT files are read inside this function.
#'
#' Workflow:
#' \enumerate{
#'   \item filter pathways by adjusted p-value and NES thresholds
#'   \item compute a Jaccard distance matrix between pathways (full gene sets)
#'   \item hierarchical clustering (default: average linkage)
#'   \item optional dendrogram branch coloring by cut height
#'   \item discrete NES heatmap (-1, 0, 1) with aligned dendrogram
#'   \item save dendrogram and heatmap as PDF files
#' }
#'
#' @param df.full A data.frame containing pathway-level GSEA results.
#'   Must include columns: `Contrast`, `Pathway`, `NES`, `p.adjust`,
#'   and `full_gene_list`.
#'   `full_gene_list` may be either a list-column of character vectors
#'   or a comma-separated character string of genes.
#' @param padj_cutoff Numeric. Adjusted p-value cutoff for pathway filtering. Default `0.05`.
#' @param threshold_positive Numeric. Positive NES threshold; values greater than this map to `+1`.
#' @param threshold_negative Numeric. Negative NES threshold; values less than this map to `-1`.
#' @param hclust_method Character. Linkage method for `stats::hclust()`. Default `"average"`.
#' @param hcut Numeric or `NA`. Height to cut dendrogram for branch coloring. Default `NA` (no coloring).
#' @param branch_palette Character vector of colors for dendrogram branches when `hcut` is specified.
#' @param branches_lwd Numeric. Line width for colored dendrogram branches.
#' @param row_split Character. `"none"` (default) or `"by_hcut"`. If `"by_hcut"` and `hcut` provided,
#'   heatmap rows are split by dendrogram clusters.
#' @param cluster_columns Logical. Whether to cluster heatmap columns (contrasts). Default `FALSE`.
#' @param heatmap_width_cm Numeric. Width of the heatmap body in centimeters.
#' @param row_dend_width_cm Numeric. Width of the row dendrogram in centimeters.
#' @param row_label_fontsize Numeric. Font size for pathway (row) labels.
#' @param row_label_bold Logical. Whether pathway labels are bold.
#' @param row_names_max_width_cm Numeric. Maximum width allowed for row labels (cm).
#' @param cell_border_col Character. Color of cell borders in the heatmap.
#' @param cell_border_lwd Numeric. Line width of cell borders.
#' @param pdf_heatmap Character. Output heatmap PDF filename.
#' @param pdf_dendrogram Character. Output dendrogram PDF filename.
#' @param pdf_width Numeric. Width of output PDFs (inches).
#' @param pdf_height Numeric. Height of output PDFs (inches).
#' @param plot_dendrogram Logical. Whether to save a standalone dendrogram PDF.
#' @param show_progress Logical. Whether to print progress during distance computation.
#'
#' @return Invisibly returns a list with components:
#' \describe{
#'   \item{df}{Filtered pathway data used for plotting}
#'   \item{dmat}{Jaccard distance matrix}
#'   \item{hc}{`hclust` object}
#'   \item{dend}{Dendrogram derived from `hc`}
#'   \item{dend_col}{Colored dendrogram (if `hcut` is used; otherwise `NULL`)}
#'   \item{nes_mat}{Numeric NES matrix (pathway x contrast)}
#'   \item{nes_disc}{Discrete NES matrix with values -1, 0, 1}
#'   \item{heatmap}{`ComplexHeatmap` heatmap object}
#' }
#'
#' @details
#' Jaccard distance is defined as:
#' \deqn{1 - \frac{|A \cap B|}{|A \cup B|}}
#'
#' Discrete NES mapping:
#' \itemize{
#'   \item `NES > threshold_positive -> +1`
#'   \item `NES < threshold_negative -> -1`
#'   \item otherwise `-> 0`
#' }
#'
#' @seealso
#' `ComplexHeatmap::Heatmap()`, `dendextend::color_branches()`, `stats::hclust()`
#'
#' @importFrom stats hclust as.dist cutree
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @importFrom grid unit gpar
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr filter select distinct mutate group_by summarise count n_distinct
#' @importFrom tibble column_to_rownames
#' @importFrom stats as.dendrogram
#' @importFrom utils capture.output flush.console head
#'
#' @export
make_gsea_jaccard_heatmap <- function(
    df.full,
    padj_cutoff = 0.05,
    threshold_positive = 2,
    threshold_negative = -2,
    hclust_method = "average",
    hcut = NA_real_,
    branch_palette = NULL,
    branches_lwd = 1.2,
    row_split = c("none", "by_hcut"),
    cluster_columns = FALSE,
    heatmap_width_cm = 2.5,
    row_dend_width_cm = 12,
    row_label_fontsize = 6,
    row_label_bold = TRUE,
    row_names_max_width_cm = 20,
    cell_border_col = "black",
    cell_border_lwd = 0.4,
    pdf_heatmap = "NES_heatmap.pdf",
    pdf_dendrogram = "full_gene_list_dendrogram.pdf",
    pdf_width = 14,
    pdf_height = 15,
    plot_dendrogram = TRUE,
    show_progress = FALSE
) {
  # ---- minimal dependency checks (helpful because ComplexHeatmap is Bioc) ----
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Install via BiocManager::install('ComplexHeatmap').")
  }
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("Package 'dendextend' is required. Install via install.packages('dendextend').")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Install via install.packages('circlize').")
  }

  row_split <- match.arg(row_split)

  # ---- required columns check ----
  if (!is.data.frame(df.full)) stop("df.full must be a data.frame.")
  req_cols <- c("Contrast", "Pathway", "NES", "p.adjust", "full_gene_list")
  missing_cols <- setdiff(req_cols, colnames(df.full))
  if (length(missing_cols) > 0) {
    stop("df.full is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ---- helpers ----
  split_genes <- function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return(character(0))
    genes <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
    genes <- trimws(genes)
    genes <- genes[nzchar(genes)]
    unique(genes)
  }

  gene_signature <- function(g) {
    # stable, comparable signature of a gene set
    paste(sort(unique(g)), collapse = "|")
  }

  # ---- default branch palette ----
  if (is.null(branch_palette)) {
    branch_palette <- rev(c(
      "tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
      "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
      "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"
    ))
  }

  # ---- filter ----
  df <- df.full %>%
    dplyr::filter(.data$p.adjust < padj_cutoff) %>%
    dplyr::filter(.data$NES > threshold_positive | .data$NES < threshold_negative)

  if (nrow(df) == 0) {
    stop("After filtering (padj_cutoff + NES thresholds), df has 0 rows. Try relaxing thresholds.")
  }

  # ---- ensure full_gene_list is list-column of character vectors ----
  if (!is.list(df$full_gene_list)) {
    df <- df %>% dplyr::mutate(full_gene_list = lapply(.data$full_gene_list, split_genes))
  } else {
    ok <- vapply(df$full_gene_list, is.character, logical(1))
    if (!all(ok)) stop("full_gene_list must be a list-column of character vectors OR a comma-separated string.")
  }

  # ---- assumption checks ----
  # 1) no duplicate (Pathway, Contrast)
  dup_pairs <- df %>% dplyr::count(.data$Pathway, .data$Contrast, name = "n") %>% dplyr::filter(.data$n > 1)
  if (nrow(dup_pairs) > 0) {
    stop("Found duplicated (Pathway, Contrast) rows. Example:\n",
         paste0(capture.output(print(head(dup_pairs, 10))), collapse = "\n"))
  }

  # 2) each Pathway has one identical gene set across contrasts
  gene_sig <- vapply(df$full_gene_list, gene_signature, character(1))

  gene_variation <- df %>%
    dplyr::mutate(.gene_sig = gene_sig) %>%
    dplyr::group_by(.data$Pathway) %>%
    dplyr::summarise(n_gene_sets = dplyr::n_distinct(.data$.gene_sig), .groups = "drop") %>%
    dplyr::filter(.data$n_gene_sets > 1)

  if (nrow(gene_variation) > 0) {
    stop("Some Pathways have >1 distinct full_gene_list across contrasts. Example:\n",
         paste0(capture.output(print(head(gene_variation, 10))), collapse = "\n"))
  }

  # ---- pathway gene sets (one per pathway) ----
  pathway_sets <- df %>%
    dplyr::mutate(.gene_sig = gene_sig) %>%
    dplyr::distinct(.data$Pathway, .data$.gene_sig, .data$full_gene_list) %>%
    dplyr::select(-.data$.gene_sig)

  gene_sets <- pathway_sets$full_gene_list
  names(gene_sets) <- pathway_sets$Pathway
  pathways <- names(gene_sets)
  n <- length(pathways)

  # ---- distance matrix ----
  dmat <- matrix(0, nrow = n, ncol = n, dimnames = list(pathways, pathways))

  if (show_progress) {
    total <- n * (n + 1) / 2
    done <- 0L
    step <- max(1L, floor(total / 100))
  }

  for (i in seq_len(n)) {
    for (j in i:n) {
      d <- jaccard_dist(gene_sets[[i]], gene_sets[[j]])
      dmat[i, j] <- d
      dmat[j, i] <- d

      if (show_progress) {
        done <- done + 1L
        if (done %% step == 0L || done == total) {
          cat(sprintf("\r%d / %d (%.1f%%) completed", done, total, 100 * done / total))
          flush.console()
        }
      }
    }
  }
  if (show_progress) cat("\n")

  hc <- stats::hclust(stats::as.dist(dmat), method = hclust_method)

  # ---- optional: save dendrogram PDF ----
  if (isTRUE(plot_dendrogram)) {
    grDevices::pdf(pdf_dendrogram, width = pdf_width, height = pdf_height, useDingbats = FALSE)
    graphics::plot(hc,
                   main = "Pathway clustering (Jaccard distance of full gene set)",
                   xlab = "", sub = "", cex = 0.6)
    grDevices::dev.off()
  }

  # ---- optional: color branches by cut height ----
  dend <- as.dendrogram(hc)
  k <- NA_integer_
  dend_col <- NULL

  if (!is.na(hcut)) {
    k <- length(unique(stats::cutree(hc, h = hcut)))
    if (k > length(branch_palette)) {
      stop("Need at least ", k, " colors in branch_palette, but got ", length(branch_palette), ".")
    }
    dend_col <- dendextend::color_branches(dend, k = k, col = branch_palette[seq_len(k)])
    dend_col <- dendextend::set(dend_col, "branches_lwd", branches_lwd)
    row_clust <- dend_col
  } else {
    row_clust <- dend
  }

  # ---- NES matrix (rows=Pathway, cols=Contrast), missing => 0 ----
  nes_mat <- df %>%
    dplyr::select(.data$Pathway, .data$Contrast, .data$NES) %>%
    tidyr::pivot_wider(names_from = .data$Contrast, values_from = .data$NES, values_fill = 0) %>%
    tibble::column_to_rownames("Pathway") %>%
    as.matrix()

  # Align matrix rows to distance-matrix order (critical)
  stopifnot(setequal(rownames(dmat), rownames(nes_mat)))
  nes_mat <- nes_mat[rownames(dmat), , drop = FALSE]
  stopifnot(identical(rownames(nes_mat), rownames(dmat)))

  # ---- discrete NES mapping ----
  nes_disc <- matrix(0, nrow = nrow(nes_mat), ncol = ncol(nes_mat), dimnames = dimnames(nes_mat))
  nes_disc[nes_mat > threshold_positive] <- 1
  nes_disc[nes_mat < threshold_negative] <- -1

  col_fun_disc <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  row_gp <- if (isTRUE(row_label_bold)) {
    grid::gpar(fontsize = row_label_fontsize, fontface = "bold")
  } else {
    grid::gpar(fontsize = row_label_fontsize)
  }

  # ---- heatmap ----
  hm <- ComplexHeatmap::Heatmap(
    nes_disc,
    name = "NES",
    col = col_fun_disc,

    rect_gp = grid::gpar(col = cell_border_col, lwd = cell_border_lwd),

    row_names_gp = row_gp,
    row_names_max_width = grid::unit(row_names_max_width_cm, "cm"),

    cluster_rows = row_clust,
    row_dend_reorder = FALSE,
    show_row_dend = TRUE,
    row_dend_width = grid::unit(row_dend_width_cm, "cm"),

    cluster_columns = isTRUE(cluster_columns),
    width = grid::unit(heatmap_width_cm, "cm"),

    row_names_side = "left",
    column_names_side = "top",
    column_names_gp = grid::gpar(fontsize = 10),

    row_split = if (row_split == "by_hcut" && !is.na(hcut)) k else NULL,

    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      labels = c(
        paste0("< ", threshold_negative),
        paste0("[", threshold_negative, ", ", threshold_positive, "]"),
        paste0("> ", threshold_positive)
      )
    )
  )

  # ---- save heatmap PDF ----
  grDevices::pdf(pdf_heatmap, width = pdf_width, height = pdf_height, useDingbats = FALSE)
  ComplexHeatmap::draw(hm)
  grDevices::dev.off()

  invisible(list(
    df = df,
    dmat = dmat,
    hc = hc,
    dend = dend,
    dend_col = dend_col,
    nes_mat = nes_mat,
    nes_disc = nes_disc,
    heatmap = hm
  ))
}
