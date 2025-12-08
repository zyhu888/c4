#' Covariate-Assisted Spectral Clustering (CASC)
#'
#' Perform covariate-assisted spectral clustering (CASC) that integrates
#' network connectivity and node covariates. The method searches over a
#' range of tuning parameters (\eqn{\alpha}) to balance structural and
#' covariate information and selects the best clustering result.
#'
#' @importFrom stats kmeans
#' @importFrom MASS mvrnorm
#'
#' @param Adj An adjacency matrix of the network (symmetric, non-negative).
#' @param Covariate A node covariate matrix where rows correspond to nodes
#' and columns to covariates.
#' @param K The number of clusters to partition the nodes into (integer).
#' @param n Number of \eqn{\alpha} values to evaluate in the grid search range. Default is 5.
#' @param ... Additional arguments passed to \code{kmeans()}, such as
#' \code{nstart}, \code{iter.max}, or initialization settings.
#'
#' @details
#' Covariate-assisted spectral clustering (CASC) is a community detection
#' algorithm for networks with node covariates, proposed by
#' Binkiewicz, Vogelstein, and Rohe (2017).
#' Let \eqn{X \in \mathbb{R}^{n \times p}} be the covariate matrix,
#' where each row corresponds to a node’s covariate vector.
#' Define the regularized graph Laplacian as
#' \deqn{L_\tau = D_\tau^{-1/2} A D_\tau^{-1/2},}
#' where \eqn{D_\tau = D + \tau I} and \eqn{D} is the diagonal degree matrix of Adj matrix \eqn{A}.
#' CASC then forms a balanced matrix:
#' \deqn{\tilde{L}(\alpha) = L_\tau L_\tau + \alpha X X^\top,}
#' where \eqn{\alpha} is a tuning parameter controlling the relative weight
#' of covariate information.
#' Clustering is performed by applying k-means on the first \eqn{K}
#' leading eigenvectors of \eqn{\tilde{L}(\alpha)}.
#'
#' @return An object of class \code{"CASC"} containing:
#' \item{alpha}{The selected tuning parameter value.}
#' \item{alpha range}{Lower and Upper bound of the search range for \eqn{\alpha}.}
#' \item{cluster}{Cluster membership assignments for each node.}
#'
#' @references
#' Binkiewicz, N., Vogelstein, J. T., & Rohe, K. (2017).
#' Covariate-assisted spectral clustering. *Biometrika*, 104(2), 361–377.
#'
#' @examples
#' library(MASS)
#' set.seed(123)
#'
#' n <- 200; k <- 4; d <- 10; se <- 2
#' z <- rep(1:k, each = n/k)
#'
#' # --- SBM adjacency ---
#' p <- matrix(0.2, k, k); diag(p) <- 0.3
#' A <- (matrix(runif(n^2), n, n) < p[z, z]) * 1
#' A[lower.tri(A)] <- t(A)[lower.tri(A)]
#' diag(A) <- 0
#'
#' # --- Covariates ---
#' means <- matrix(c(
#'   rep(d / sqrt(8), 3),
#'   rep(c(-d, -d, d) / sqrt(8), 1),
#'   rep(c(-d, d, -d) / sqrt(8), 1),
#'   rep(c(d, -d, -d) / sqrt(8), 1)
#' ), nrow = k, byrow = TRUE)
#' X <- do.call(rbind, lapply(1:k, \(i) MASS::mvrnorm(n/k, means[i, ], diag(se^2, 3))))
#'
#' # Run CASC clustering
#' result_CASC <- CASC(A, X, K = 4, iter.max = 100, nstart = 10)
#' result_CASC
#'
#' @export
CASC <- function(Adj, Covariate, K, n = 5, ...) {

  # Normalized Laplacian
  D <- diag(rowSums(Adj) + mean(rowSums(Adj)))
  D_inv_sqrt <- diag(diag(D)^(-0.5))
  V <- D_inv_sqrt %*% Adj %*% D_inv_sqrt

  # Eigen decomposition
  net_eigen <- eigen(V %*% V, symmetric = TRUE)
  ca <- Covariate %*% t(Covariate)
  ca_eigen <- eigen(ca, symmetric = TRUE)

  # candidate alpha
  R <- ncol(Covariate)
  alpha_upper <- ifelse(R <= K,
                        net_eigen$values[1] / ca_eigen$values[R],
                        net_eigen$values[1] / (ca_eigen$values[K] - ca_eigen$values[K + 1]))
  alpha_lower <- (net_eigen$values[K] - net_eigen$values[K + 1]) / ca_eigen$values[1]
  alphagrid <- seq(alpha_lower, alpha_upper, length.out = n)

  # cluster for each alpha
  within_ss <- numeric(n)
  clusters <- matrix(0, n, nrow(Adj))

  for (i in 1:n) {
    eigen_decomp <- eigen(V %*% V + alphagrid[i] * ca, symmetric = TRUE)
    X <- eigen_decomp$vectors[, 1:K]
    row_norms <- sqrt(rowSums(X^2))
    valid_idx <- row_norms > 0
    X[valid_idx, ] <- X[valid_idx, ] / row_norms[valid_idx]
    kmeans_result <- kmeans(X[valid_idx, ], K,...)
    clusters[i, valid_idx] <- kmeans_result$cluster
    within_ss[i] <- kmeans_result$tot.withinss
  }
  # return result
  best_idx <- which.min(within_ss)
  result <- list(
    cluster = clusters[best_idx, ],
    K = K,
    alpha = alphagrid[best_idx],
    alphalower = alpha_lower,
    alphaupper = alpha_upper
  )
  class(result) <- "CASC"
  return(result)
}

