print_cluster_result <- function(x, title, ...) {
  cat(title, "\n")
  cat("------------------------------------------------\n")

  ## Selected alpha
  cat("Selected alpha: ", x$alpha, "\n")

  ## Alpha range
  cat("Alpha range:  [", x$alphalower, ", ", x$alphaupper, "]\n", sep = "")

  ## Number of clusters (K)
  cat("Number of clusters (K): ", x$K, "\n")

  ## Cluster assignments
  cat("Cluster assignments:\n")
  print(x$cluster)

  invisible(x)
}

#' @export
print.CASC <- function(x, ...) {
  print_cluster_result(
    x,
    title = "Covariate-Assisted Spectral Clustering (CASC)",
    ...
  )
}

#' @export
print.C4 <- function(x, ...) {
  print_cluster_result(
    x,
    title = "Covariate Connectivity Combined Clustering (C4)",
    ...
  )
}
