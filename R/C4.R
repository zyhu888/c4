#' C4: Covariate Connectivity Combined Clustering
#'
#' Perform C4 clustering, an adaptive spectral clustering method that integrates
#' network connectivity with node covariates for unknown number of communities.
#' The algorithm searches over a sequence of tuning parameters (\eqn{\alpha}) to
#' balance adjacency and covariate similarity, and selects the best result based
#' on silhouette score.
#'
#' @importFrom stats dist as.dist kmeans
#' @importFrom cluster silhouette
#' @importFrom MASS mvrnorm
#'
#' @param adj_matrix Adjacency matrix of the network (symmetric, non-negative, with zero diagonal).
#' @param sim_matrix Similarity matrix derived from covariates (symmetric, non-negative, with zero diagonal).
#' @param K Number of clusters, or a range of integers. If a range is given,
#' the eigengap heuristic is used to select the best \code{K}. Default is \code{2:(n-1)}.
#' @param alpha_seq A numeric vector of candidate alpha values. Default is \code{seq(0, 1, by = 0.1)}.
#' @param epsilon A small positive constant added to avoid division by zero
#' when computing inverse distances. Default is \code{1e-6}.
#'
#' @details
#' The method forms a convex combination of adjacency and covariate similarity matrices:
#' \deqn{(1-\alpha)W + \alpha S,}
#' where \eqn{W} is the adjacency matrix and \eqn{S} is the covariate similarity.
#' For each \eqn{\alpha}, the normalized Laplacian is computed, eigen-decomposition
#' is performed, and k-means is applied on the spectral embedding. The solution
#' with the highest silhouette score is returned.
#'
#' @return A list containing:
#' \item{alpha}{The selected alpha value.}
#' \item{K}{The selected number of clusters.}
#' \item{cluster}{Cluster membership assignments for each node.}
#'
#' @references Work in progress (manuscript in preparation).#'
#'
#' @examples
#' library(MASS)
#' set.seed(123)
#'
#' # Generate covariates (200 nodes, 4 communities, 3 features per node)
#' # Each community has a different mean vector in 3D space
#' n <- 200; k <- 4; d <- 10; se <- 3
#' means <- matrix(c(
#'   d / sqrt(8),  d / sqrt(8),  d / sqrt(8),
#'   -d / sqrt(8), -d / sqrt(8),  d / sqrt(8),
#'   -d / sqrt(8),  d / sqrt(8), -d / sqrt(8),
#'   d / sqrt(8), -d / sqrt(8), -d / sqrt(8)
#' ), nrow = 4, byrow = TRUE)
#' cov_matrix <- diag(se^2, 3)
#' community_sizes <- rep(n / k, k)
#'
#' # Simulate multivariate normal features for each community
#' X <- do.call(rbind, lapply(1:k, function(c) {
#'   mvrnorm(community_sizes[c], means[c, ], cov_matrix)
#' }))
#'
#' # Compute similarity matrix from covariates (Euclidean → 1/distance)
#' D <- as.matrix(dist(X))
#' S <- 1 / D
#' diag(S) <- 0
#'
#' # Generate adjacency matrix from a weighted SBM
#' W <- gen_weighted_sbm(
#'   node_num = 200,
#'   cluster_size = rep(50, 4),
#'   win_cluster_den = 0.3,
#'   win_cluster_dist = "runif",
#'   win_cluster_par = list(min = 1, max = 2),
#'   btw_cluster_den = 0.2,
#'   btw_cluster_dist = "runif",
#'   btw_cluster_par = list(min = 1, max = 2)
#' )$adjacency_matrix
#'
#' # Run C4 clustering
#' result_C4 <- C4(W, S, alpha_seq = seq(0, 1, by = 0.1))
#' result_C4
#'
#' @export
C4 <- function(adj_matrix, sim_matrix, K = NULL, alpha_seq = seq(0, 1, by = 0.1),
               epsilon = 1e-6) {
  n <- nrow(adj_matrix)

  # If K is not provided, set default K = 2:n-1
  if (is.null(K)) {
    K <- 2:(n-1)
  }

  # Initialize best metrics
  best_silhouette <- -Inf
  best_alpha <- NA
  best_K <- NA
  best_cluster <- NULL

  # Compute scaling factor based on the sums of the adjacency and similarity matrices
  diag(sim_matrix) <- 0
  scale_factor <- sum(adj_matrix) / sum(sim_matrix)
  sim_matrix_scaled <- sim_matrix * scale_factor

  # Loop over all candidate alpha values
  for (alpha in alpha_seq) {
    # Combined matrix
    combined_matrix <- (1 - alpha) * adj_matrix + alpha * sim_matrix_scaled
    combined_matrix <- as.matrix(combined_matrix)

    # Normalized Laplacian
    D <- diag(rowSums(combined_matrix))
    D_inv_sqrt <- solve(sqrt(D))
    V <- D_inv_sqrt %*% combined_matrix %*% D_inv_sqrt
    L <- diag(nrow(D)) - V

    # Eigen-decomposition
    eigen_decomp <- eigen(L, symmetric = TRUE)
    eigenvalues <- sort(eigen_decomp$values)
    eigenvectors <- eigen_decomp$vectors[, order(eigen_decomp$values)]

    # Eigengap heuristic if K is a range
    if (length(K) > 1) {
      eigengap <- diff(eigenvalues)
      eigengap_subset <- eigengap[(min(K)):(max(K))]
      best_k <- which.max(eigengap_subset) + min(K) - 1
    } else {
      best_k <- K
    }

    # Spectral embedding
    X <- eigenvectors[, 1:best_k]
    Y <- X / sqrt(rowSums(X^2))

    # k-means clustering
    kmeans_result <- kmeans(Y, centers = best_k, nstart = 10)
    membership <- kmeans_result$cluster

    # Silhouette score on inverse-distance graph
    dist_graph <- as.dist(1 / (combined_matrix + epsilon))
    sil <- cluster::silhouette(membership, dist_graph)
    avg_sil <- mean(sil[, "sil_width"])

    # Update best result if silhouette improves
    if (avg_sil > best_silhouette) {
      best_silhouette <- avg_sil
      best_alpha <- alpha
      best_K <- best_k
      best_cluster <- membership
    }
  }

  return(list(
    alpha = best_alpha,
    K = best_K,
    cluster = best_cluster
  ))
}
