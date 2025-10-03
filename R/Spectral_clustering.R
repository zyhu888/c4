#' Spectral Clustering with Eigengap Heuristic
#'
#' Perform traditional spectral clustering on an adjacency matrix. If the number of clusters \code{K}
#' is not specified, the function applies the eigengap heuristic (Ng, Jordan, & Weiss, 2002)
#' to determine the optimal number of communities.
#'
#'@importFrom stats kmeans
#'
#' @param adj_matrix A square adjacency matrix representing the network (numeric, symmetric, non-negative).
#' @param K An integer or vector of integers specifying the number of clusters.
#' If \code{NULL}, defaults to \code{2:(n-1)}, where \code{n} is the number of nodes.
#' If a range is provided, the eigengap heuristic is used to select the best \code{K}.
#'
#' @details
#' The algorithm computes the normalized graph Laplacian, performs eigen-decomposition,
#' and applies k-means clustering on the top \code{K} eigenvectors after row normalization.
#'
#' @return A list containing:
#' \item{K}{The selected number of clusters.}
#' \item{kmeans_result}{The \code{kmeans} result object.}
#' \item{membership}{A vector of cluster assignments for each node.}
#'
#' @references
#' Ng, A. Y., Jordan, M. I., & Weiss, Y. (2002).
#' On spectral clustering: Analysis and an algorithm.
#' *Advances in Neural Information Processing Systems (NIPS)*.
#'
#' von Luxburg, U. (2007).
#' A tutorial on spectral clustering.
#' *Statistics and Computing*, 17(4), 395–416.
#'
#' @examples
#' # Example: simple adjacency matrix
#' adj <- matrix(c(0,1,1,0,
#'                 1,0,0,1,
#'                 1,0,0,1,
#'                 0,1,1,0), nrow = 4, byrow = TRUE)
#'
#' result <- Spectral_clustering(adj, K = 2)
#' result$membership
#'
#' @export
Spectral_clustering <- function(adj_matrix, K = NULL) {
  n <- nrow(adj_matrix)

  # If K is not specified, default to 2:n-1
  if (is.null(K)) {
    K <- 2:(n-1)
  }
  adj_matrix <- as.matrix(adj_matrix)

  # Compute degree matrix and normalized matrix
  D <- diag(rowSums(adj_matrix))
  D_inv_sqrt <- solve(sqrt(D))
  V <- D_inv_sqrt %*% adj_matrix %*% D_inv_sqrt
  L <- diag(nrow(D)) - V

  # Perform eigen-decomposition of the normalized Laplacian matrix
  eigen_decomp <- eigen(L, symmetric = TRUE)
  eigenvalues <- sort(eigen_decomp$values)
  eigenvectors <- eigen_decomp$vectors[, order(eigen_decomp$values)]

  # If K is a range, use the eigengap heuristic to select the best K
  if (length(K) > 1) {
    eigengap <- diff(eigenvalues)
    eigengap_subset <- eigengap[(min(K)):(max(K))]
    best_k <- which.max(eigengap_subset) + min(K) - 1
  } else {
    best_k <- K
  }

  # Select the first best_k eigenvectors
  X <- eigenvectors[,1:best_k]
  Y <- X / sqrt(rowSums(X^2))

  # Apply k-means clustering on the normalized eigenvector matrix
  kmeans_result <- kmeans(Y, centers = best_k, nstart = 10)
  membership <- kmeans_result$cluster

  return(list(
    K = best_k,
    kmeans_result = kmeans_result,
    membership = kmeans_result$cluster
  ))
}
