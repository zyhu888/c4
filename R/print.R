#' @export
print.CASC <- function(x, ...) {
  cat("Covariate-Assisted Spectral Clustering (CASC)\n")
  cat("------------------------------------------------\n")
  cat("Selected alpha: ", x$alpha, "\n")
  cat("Alpha range:  [", x$alphalower, ", ", x$alphaupper, "]\n", sep = "")
  cat("Cluster assignments:\n")
  print(x$cluster)

  invisible(x)
}

#' @export
print.C4 <- function(x, ...) {
  cat("Covariate Connectivity Combined Clustering (C4)\n")
  cat("------------------------------------------------\n")
  cat("Selected alpha: ", x$alpha, "\n")
  cat("Number of clusters (K): ", x$K, "\n")
  cat("Cluster assignments:\n")
  print(x$cluster)

  invisible(x)
}
