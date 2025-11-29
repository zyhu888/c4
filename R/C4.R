#' Covariate Connectivity Combined Clustering (C⁴)
#'
#' Perform Covariate Connectivity Combined Clustering (C⁴), an adaptive spectral
#' clustering method that integrates network connectivity with node covariates
#' for unknown number of communities. The algorithm searches over a sequence of
#' tuning parameters (\eqn{\alpha}) to balance adjacency and covariate similarity,
#' and selects the best result based on silhouette score.
#'
#' @importFrom stats dist as.dist kmeans
#' @importFrom cluster silhouette
#' @importFrom MASS mvrnorm
#'
#' @param adj Adjacency matrix of the network (symmetric, non-negative, with zero diagonal).
#' @param sim Optional similarity matrix derived from covariates (symmetric, non-negative, with zero diagonal).
#' If not provided, the function sets \eqn{\alpha = 0} to perform standard spectral clustering
#' based solely on Adjacency matrix.
#' @param K Number of clusters, or a range of integers. If a range is given,
#' the eigengap heuristic is used to select the best \code{K}. Default is \code{2:(n-1)}.
#' @param alphagrid A numeric vector of candidate alpha values. Default is \code{seq(0, 1, by = 0.1)}. Ignored when \code{sim = NULL}.
#' @param ... Additional arguments passed to \code{kmeans()}, such as
#' \code{nstart}, \code{iter.max}, or initialization settings.
#'
#' @details
#' The method forms a combination of adjacency and covariate similarity matrices:
#' \deqn{(1-\alpha)W + \alpha S,}
#' where \eqn{W} is the adjacency matrix and \eqn{S} is the covariate similarity.
#' For each \eqn{\alpha}, the normalized Laplacian is computed, eigen-decomposition
#' is performed, number of clusters is determined by eigengap, and k-means is applied
#' on the spectral embedding. The clustering result with the highest silhouette score is returned.
#'
#' @return An object of class \code{"C4"} containing:
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
#' n <- 200; k <- 4; d <- 10; se <- 2
#' z <- rep(1:k, each = n/k)
#'
#' # --- SBM adjacency ---
#' p <- matrix(0.2, k, k); diag(p) <- 0.3
#' A <- (matrix(runif(n^2), n, n) < p[z, z]) * 1
#' A[lower.tri(A)] <- t(A)[lower.tri(A)]
#' diag(A) <- 0
#' # --- Covariates ---
#' means <- matrix(c(
#'   rep(d / sqrt(8), 3),
#'   rep(c(-d, -d, d) / sqrt(8), 1),
#'   rep(c(-d, d, -d) / sqrt(8), 1),
#'   rep(c(d, -d, -d) / sqrt(8), 1)
#' ), nrow = k, byrow = TRUE)
#' X <- do.call(rbind, lapply(1:k, \(i) MASS::mvrnorm(n/k, means[i, ], diag(se^2, 3))))
#'
#' # Compute similarity matrix from covariates (Euclidean → 1/distance)
#' D <- as.matrix(dist(X))
#' S <- 1 / D
#' diag(S) <- 0
#'
#' # C4 clustering on both adj and sim
#' result_C4 <- C4(A, S, iter.max = 100, nstart = 10)
#' # spectrum clustering on adj only
#' result_spe <- C4(A)
#'
#' result_C4
#' result_spe
#'
#' @export
C4 <- function(adj, sim = NULL, K = NULL, alphagrid = seq(0, 1, by = 0.1),...) {
  # Checks for adjacency matrix
  adj <- as.matrix(adj)

  if (nrow(adj) != ncol(adj)) {
    stop("Adjacency matrix must be square!", call. = FALSE)
  }

  if (any(adj < 0, na.rm = TRUE)) {
    stop("Adjacency matrix must be non-negative!", call. = FALSE)
  }

  if (!isTRUE(all.equal(adj, t(adj)))) {
    stop("Adjacency matrix must be symmetric.", call. = FALSE)
  }

  if (any(diag(adj) != 0, na.rm = TRUE)) {
    warning("Adjacency matrix has non-zero diagonal; setting diagonal entries to 0!",
            call. = FALSE)
    diag(adj) <- 0
  }

  n <- nrow(adj)

  if (is.null(K)) {
    K <- 2:(n - 1)
  }

  # If sim is NULL, force alpha = 0 (standard spectral clustering)
  if (is.null(sim)) {
    alphagrid <- 0
    message("No similarity matrix is provided; running standard spectral clustering (alpha = 0)!")
  } else {
    ## Checks for sim
    sim <- as.matrix(sim)

    if (nrow(sim) != ncol(sim) || nrow(sim) != n) {
      stop("Similarity matrix must be square and have the same dimension as adjacency matrix",
           call. = FALSE)
    }

    if (!isTRUE(all.equal(sim, t(sim)))) {
      stop("Similarity matrix must be symmetric!", call. = FALSE)
    }

    if (any(diag(sim) != 0, na.rm = TRUE)) {
      warning("Similarity matrix has non-zero diagonal; setting diagonal entries to 0.",
              call. = FALSE)
      diag(sim) <- 0
    }

    ## Compute scaling factor only if sim is provided
    if (sum(sim) == 0) {
      stop("Similarity matrix 'sim' sums to zero; cannot compute scaling factor.",
           call. = FALSE)
    }
    scale_factor <- sum(adj) / sum(sim)
    sim_scaled   <- sim * scale_factor
  }

  # Initialize best metrics
  best_silhouette <- -Inf
  best_alpha <- NA
  best_K <- NA
  best_cluster <- NULL

  # Loop over all candidate alpha values
  for (alpha in alphagrid) {
    if (is.null(sim)) {
      # Pure spectral clustering: use only adjacency matrix
      combined_matrix <- adj
    } else {
      # Combined matrix
      combined_matrix <- (1 - alpha) * adj + alpha * sim_scaled
    }
    combined_matrix <- as.matrix(combined_matrix)

    # Normalized Laplacian
    D <- diag(rowSums(combined_matrix))
    D_inv_sqrt <- solve(sqrt(D))
    V <- D_inv_sqrt %*% combined_matrix %*% D_inv_sqrt
    L <- diag(nrow(D)) - V

    # Eigen decomposition
    eigen_decomp <- eigen(L, symmetric = TRUE)
    eigenvalues <- sort(eigen_decomp$values)
    eigenvectors <- eigen_decomp$vectors[, order(eigen_decomp$values)]

    # Eigengap heuristic
    if (length(K) > 1) {
      eigengap <- diff(eigenvalues)
      eigengap_subset <- eigengap[(min(K)):(max(K))]
      best_k <- which.max(eigengap_subset) + min(K) - 1
    } else {
      best_k <- K
    }

    # k-means clustering
    X <- eigenvectors[, 1:best_k]
    Y <- X / sqrt(rowSums(X^2))
    kmeans_result <- kmeans(Y, centers = best_k, ...)
    membership <- kmeans_result$cluster

    # Only compute silhouette if sim is provided (for alpha selection)
    if (!is.null(sim)) {
      dist_graph <- as.dist(1 / (combined_matrix + 1e-6))
      sil <- cluster::silhouette(membership, dist_graph)
      avg_sil <- mean(sil[, "sil_width"])

      if (avg_sil > best_silhouette) {
        best_silhouette <- avg_sil
        best_alpha <- alpha
        best_K <- best_k
        best_cluster <- membership
      }
    } else {
      # No alpha selection needed
      best_alpha <- alpha
      best_K <- best_k
      best_cluster <- membership
      best_silhouette <- NA
    }
  }
 result <- list(
   alpha = best_alpha,
   K = best_K,
   cluster = best_cluster
 )
 class(result) <- "C4"
 return(result)
}
