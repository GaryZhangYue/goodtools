#' Compute Jaccard distance between two gene sets
#'
#' @param A Character vector of genes.
#' @param B Character vector of genes.
#' @return Numeric scalar distance between 0 and 1. A higher distance value indicates a bigger difference.
#' @export
jaccard_dist <- function(A, B) {
  A <- unique(A); B <- unique(B)
  if (length(A) == 0 && length(B) == 0) return(0)
  u <- union(A, B)
  if (length(u) == 0) return(0)
  1 - (length(intersect(A, B)) / length(u))
}
