#' Plot DEG barplot with pathway-membership dot-heatmap (GSEA leading-edge)
#'
#' For a given Contrast × cluster, this function:
#' \enumerate{
#'   \item counts how many pathways each gene appears in (from GSEA leadingEdge),
#'   \item merges in DEG log fold-change (lfc) and a precomputed padj category (padj_cat),
#'   \item draws a faceted barplot (facet = n_pathways, ordered high → low),
#'   \item draws a pathway-membership "dot heatmap" aligned underneath (binary present/absent).
#' }
#'
#' The heatmap panel is rendered as a grid of white tiles (for alignment) with colored circles
#' indicating membership in each pathway. The bar/heatmap facets share the same ordering.
#'
#' Users can further customize plots by providing ggplot components via:
#' \itemize{
#'   \item \code{...}: applied to BOTH \code{p_bar} and \code{p_heat}
#'   \item \code{bar_add}: applied only to \code{p_bar}
#'   \item \code{heat_add}: applied only to \code{p_heat}
#' }
#'
#' @param df_gsea Unnested GSEA results table. Must include columns for contrast, cluster,
#'   pathway, and a list-column of leading-edge genes (e.g. \code{leadingEdge}).
#' @param df_deg DEG table. Must include gene, contrast, log fold-change, and \code{padj_cat}.
#'   \code{padj_cat} should be computed outside this function.
#' @param contrast Character scalar. Which contrast to plot.
#' @param cluster Character scalar. Which cluster to plot (matched as character).
#'
#' @param col_contrast Column name for contrast in both tables (default \code{"Contrast"}).
#' @param col_cluster Column name for cluster in \code{df_gsea} (default \code{"Cluster"}).
#' @param col_pathway Column name for pathway (default \code{"Pathway"}).
#' @param col_leading Column name for leading-edge list-column (default \code{"leadingEdge"}).
#' @param col_gene Column name for gene (default \code{"Gene"}).
#' @param col_lfc Column name for log fold-change in \code{df_deg} (default \code{"lfc"}).
#' @param col_padj_cat Column name for padj category in \code{df_deg} (default \code{"padj_cat"}).
#'
#' @param legend_position Legend position for barplot (default \code{"top"}).
#' @param legend_direction Legend direction for barplot (default \code{"horizontal"}).
#'
#' @param heat_colors Named character vector for membership colors with names \code{"0"} and \code{"1"}
#'   (default \code{c(`0`="white", `1`="lightblue")}).
#' @param circle_size Numeric; size of membership circles in the heatmap (default 2.5).
#' @param heights Numeric length-2 vector for patchwork layout heights (default \code{c(1,2)}).
#'
#' @param x_text_angle Angle for barplot x labels (default 90).
#' @param x_text_size Size for barplot x labels (default 7).
#' @param pathway_text_size Size for pathway labels (right side) (default 8).
#'
#' @param ... Additional ggplot components (e.g. \code{theme()}, \code{scale_*()}, \code{guides()})
#'   that will be added to BOTH plots.
#' @param bar_add A list of ggplot components to add only to the barplot (default \code{list()}).
#' @param heat_add A list of ggplot components to add only to the heatmap (default \code{list()}).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{plot}: patchwork combined plot
#'   \item \code{p_bar}: barplot ggplot
#'   \item \code{p_heat}: heatmap ggplot
#'   \item \code{bar_df}: data used for barplot
#'   \item \code{heat_df}: data used for heatmap
#'   \item \code{gene_freq_deg}: merged gene frequency + DEG stats
#' }
#'
#' @examples
#' \dontrun{
#' # df_deg must already contain padj_cat
#' res <- plot_deg_gsea_bar_heat(
#'   df_gsea   = df.tmp.unnested,
#'   df_deg    = df.deg2,
#'   contrast  = "Classical Monocytes",
#'   cluster   = "1",
#'   circle_size = 3,
#'   bar_add  = list(ggplot2::theme(legend.position = "bottom")),
#'   heat_add = list(ggplot2::theme(axis.text.y.right = ggplot2::element_text(size = 10)))
#' )
#' res$plot
#' }
#'
#' @importFrom dplyr ungroup mutate filter transmute distinct count left_join arrange group_by ungroup select rename
#' @importFrom tidyr unnest replace_na
#' @importFrom ggplot2 ggplot aes geom_col geom_tile geom_point facet_grid scale_color_manual scale_y_discrete labs theme theme_bw theme_minimal element_text element_blank
#' @importFrom patchwork plot_layout
#' @importFrom rlang abort
#' @export
plot_deg_gsea_bar_heat <- function(
    df_gsea,
    df_deg,
    contrast,
    cluster,
    col_contrast = "Contrast",
    col_cluster  = "Cluster",
    col_pathway  = "Pathway",
    col_leading  = "leadingEdge",
    col_gene     = "Gene",
    col_lfc      = "lfc",
    col_padj_cat = "padj_cat",
    legend_position  = "top",
    legend_direction = "horizontal",
    heat_colors = c(`0` = "white", `1` = "lightblue"),
    circle_size = 2.5,
    heights = c(1, 2),
    x_text_angle = 90,
    x_text_size = 7,
    pathway_text_size = 8,
    ...,
    bar_add = list(),
    heat_add = list()
) {
  # ---- deps
  for (pkg in c("dplyr", "tidyr", "ggplot2", "patchwork", "rlang")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
    }
  }

  # ---- basic input checks
  if (!is.data.frame(df_gsea)) rlang::abort("df_gsea must be a data.frame / tibble.")
  if (!is.data.frame(df_deg))  rlang::abort("df_deg must be a data.frame / tibble.")
  if (!(is.character(contrast) && length(contrast) == 1 && nzchar(contrast))) {
    rlang::abort("contrast must be a non-empty character scalar.")
  }
  if (!(is.character(cluster) && length(cluster) == 1 && nzchar(cluster))) {
    rlang::abort("cluster must be a non-empty character scalar (pass '1', '2', etc).")
  }
  if (!(is.numeric(heights) && length(heights) == 2 && all(is.finite(heights)) && all(heights > 0))) {
    rlang::abort("heights must be a numeric length-2 vector of positive values, e.g. c(1, 2).")
  }
  if (!(is.numeric(circle_size) && length(circle_size) == 1 && is.finite(circle_size) && circle_size > 0)) {
    rlang::abort("circle_size must be a positive numeric scalar.")
  }
  if (!all(c("0","1") %in% names(heat_colors))) {
    rlang::abort("heat_colors must be a named vector with names '0' and '1'.")
  }
  if (!is.list(bar_add) || !is.list(heat_add)) {
    rlang::abort("bar_add and heat_add must be lists of ggplot components, e.g. list(theme(...), scale_*()).")
  }

  # ---- column presence checks
  must_have_gsea <- c(col_contrast, col_cluster, col_pathway, col_leading)
  missing_gsea <- setdiff(must_have_gsea, names(df_gsea))
  if (length(missing_gsea) > 0) {
    rlang::abort(paste0("df_gsea is missing required columns: ", paste(missing_gsea, collapse = ", ")))
  }

  must_have_deg <- c(col_gene, col_contrast, col_lfc, col_padj_cat)
  missing_deg <- setdiff(must_have_deg, names(df_deg))
  if (length(missing_deg) > 0) {
    rlang::abort(paste0("df_deg is missing required columns: ", paste(missing_deg, collapse = ", ")))
  }

  # ---- sanity: leadingEdge list-column
  if (!is.list(df_gsea[[col_leading]])) {
    rlang::abort(paste0("df_gsea[['", col_leading, "']] must be a list-column (each row a vector of genes)."))
  }

  # ---- load namespaces locally (no attach)
  dplyr <- asNamespace("dplyr")
  tidyr <- asNamespace("tidyr")
  ggplot2 <- asNamespace("ggplot2")

  # ---- subset GSEA (convert cluster col to character for stable matching)
  gsea_sub <- df_gsea |>
    dplyr$ungroup() |>
    dplyr$mutate(.cluster_chr = as.character(.data[[col_cluster]])) |>
    dplyr$filter(.data[[col_contrast]] == contrast, .cluster_chr == cluster) |>
    dplyr$transmute(
      Pathway = .data[[col_pathway]],
      genes   = .data[[col_leading]]
    )

  if (nrow(gsea_sub) == 0) {
    rlang::abort(paste0(
      "No rows found in df_gsea for contrast='", contrast, "' and cluster='", cluster, "'."
    ))
  }

  # ---- gene frequency across pathways
  gene_freq <- gsea_sub |>
    tidyr$unnest(genes) |>
    dplyr$filter(!is.na(genes), genes != "") |>
    dplyr$distinct(Pathway, genes) |>
    dplyr$count(genes, name = "n_pathways") |> # count is shorthand for group_by(...) %>% summarise(n = n())
    dplyr$rename(Gene = genes)

  if (nrow(gene_freq) == 0) {
    rlang::abort("After unnesting leadingEdge genes, no genes remained (all NA/empty?).")
  }

  # ---- subset DEG + join
  deg_sub <- df_deg |>
    dplyr$ungroup() |>
    dplyr$filter(.data[[col_contrast]] == contrast) |>
    dplyr$transmute(
      Gene = as.character(.data[[col_gene]]),
      lfc_raw = .data[[col_lfc]],
      padj_cat = .data[[col_padj_cat]]
    )

  ## ---- 1. Strict numeric coercion check for lfc
  lfc_num <- suppressWarnings(as.numeric(deg_sub$lfc_raw))

  if (any(is.na(lfc_num) & !is.na(deg_sub$lfc_raw))) {
    bad_vals <- unique(deg_sub$lfc_raw[is.na(lfc_num) & !is.na(deg_sub$lfc_raw)])
    rlang::abort(
      paste0(
        "Non-numeric values found in lfc column (cannot be coerced): ",
        paste(bad_vals, collapse = ", ")
      )
    )
  }

  deg_sub <- deg_sub |>
    dplyr$mutate(lfc = lfc_num) |>
    dplyr$select(Gene, lfc, padj_cat)

  ## ---- 2. Strict duplicate gene check
  dup_genes <- deg_sub$Gene[duplicated(deg_sub$Gene)]

  if (length(dup_genes) > 0) {
    rlang::abort(
      paste0(
        "Duplicate genes found in DEG table after filtering: ",
        paste(unique(dup_genes), collapse = ", "),
        ". Each gene must appear exactly once."
      )
    )
  }

  gene_freq_deg <- gene_freq |>
    dplyr$left_join(deg_sub, by = "Gene") |>
    dplyr$mutate(
      padj_cat = as.factor(padj_cat)
    )

  # ---- bar_df (facet order high->low, within facet sort by lfc desc)
  bar_df <- gene_freq_deg |>
    dplyr$mutate(n_pathways = as.integer(n_pathways)) |>
    dplyr$group_by(n_pathways) |>
    dplyr$arrange(dplyr$desc(lfc), .by_group = TRUE) |> # .by_group = TRUE → do the sorting inside each n_pathways group
    dplyr$mutate(Gene_ord = factor(Gene, levels = unique(Gene))) |>
    dplyr$ungroup()

  facet_levels <- sort(unique(bar_df$n_pathways), decreasing = TRUE)
  bar_df <- bar_df |>
    dplyr$mutate(n_pathways_f = factor(n_pathways, levels = facet_levels))

  # ---- membership table (Pathway x Gene)
  mem <- df_gsea |>
    dplyr$ungroup() |>
    dplyr$mutate(.cluster_chr = as.character(.data[[col_cluster]])) |>
    dplyr$filter(.data[[col_contrast]] == contrast, .cluster_chr == cluster) |>
    dplyr$transmute(
      Pathway = .data[[col_pathway]],
      Gene    = .data[[col_leading]]
    ) |>
    tidyr$unnest(Gene) |>
    dplyr$filter(!is.na(Gene), Gene != "") |>
    dplyr$distinct(Pathway, Gene) |>
    dplyr$mutate(present = 1L)

  # align heatmap x to barplot genes + facet membership
  gene_key <- bar_df |>
    dplyr$select(Gene, n_pathways, n_pathways_f, Gene_ord) |>
    dplyr$distinct() |>
    dplyr$rename(Gene_ord_bar = Gene_ord)

  pathways_all <- mem |>
    dplyr$distinct(Pathway)

  heat_df <- tidyr$expand_grid(
    Pathway = pathways_all$Pathway,
    Gene    = gene_key$Gene
  ) |>
    dplyr$left_join(mem,      by = c("Pathway", "Gene")) |>
    dplyr$left_join(gene_key, by = "Gene") |>
    dplyr$mutate(
      present      = tidyr$replace_na(present, 0L),
      present      = factor(present, levels = c(0, 1)),
      n_pathways_f = factor(n_pathways_f, levels = levels(bar_df$n_pathways_f)),
      Gene_ord     = Gene_ord_bar
    ) |>
    dplyr$select(-Gene_ord_bar)

  # ---- plots
  p_bar <- ggplot2$ggplot(bar_df, ggplot2$aes(x = Gene_ord, y = lfc, fill = padj_cat)) +
    ggplot2$geom_col(width = 0.9) +
    ggplot2$facet_grid(. ~ n_pathways_f, scales = "free_x", space = "free_x") +
    ggplot2$labs(x = NULL, y = "Log Fold Change", fill = "padj") +
    ggplot2$theme_bw() +
    ggplot2$theme(
      axis.text.x = ggplot2$element_text(
        angle = x_text_angle, vjust = 0.5, hjust = 1, size = x_text_size, face = "bold"
      ),
      panel.grid.major.x = ggplot2$element_blank(),
      legend.position = legend_position,
      legend.direction = legend_direction
    )

  # Dot-heatmap: white tiles for alignment + colored circles for membership
  p_heat <- ggplot2$ggplot(heat_df, ggplot2$aes(x = Gene_ord, y = Pathway)) +
    ggplot2$geom_tile(fill = "white", color = NA) +
    ggplot2$geom_point(ggplot2$aes(color = present), shape = 16, size = circle_size) +
    ggplot2$facet_grid(. ~ n_pathways_f, scales = "free_x", space = "free_x") +
    ggplot2$scale_color_manual(values = heat_colors, guide = "none") +
    ggplot2$scale_y_discrete(position = "right") +
    ggplot2$labs(x = NULL, y = NULL) +
    ggplot2$theme_minimal() +
    ggplot2$theme(
      axis.text.x = ggplot2$element_blank(),
      axis.ticks.x = ggplot2$element_blank(),
      panel.grid = ggplot2$element_blank(),
      axis.text.y.right = ggplot2$element_text(size = pathway_text_size),
      axis.text.y.left  = ggplot2$element_blank(),
      strip.text = ggplot2$element_blank(),
      strip.background = ggplot2$element_blank()
    )

  # ---- apply user customizations
  common_add <- list(...)
  if (length(common_add) > 0) {
    p_bar  <- p_bar  + common_add
    p_heat <- p_heat + common_add
  }
  if (length(bar_add) > 0)  p_bar  <- p_bar  + bar_add
  if (length(heat_add) > 0) p_heat <- p_heat + heat_add

  combined <- p_bar / p_heat + patchwork::plot_layout(heights = heights)

  list(
    plot = combined,
    p_bar = p_bar,
    p_heat = p_heat,
    bar_df = bar_df,
    heat_df = heat_df,
    gene_freq_deg = gene_freq_deg
  )
}

