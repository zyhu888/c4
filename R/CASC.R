#' Covariate-Assisted Spectral Clustering (CASC)
#'
#' Perform covariate-assisted spectral clustering (CASC) that integrates
#' network connectivity and node covariates. The method searches over a
#' range of tuning parameters (\eqn{\alpha}) to balance structural and
#' attribute information and selects the best clustering result.
#'
#'@importFrom stats kmeans
#'
#' @param Adj An adjacency matrix of the network (symmetric, non-negative).
#' @param Covariate A node covariate matrix where rows correspond to nodes
#' and columns to covariates.
#' @param K The number of clusters to partition the nodes into (integer).
#' @param alphan Number of candidate \eqn{\alpha} values to evaluate. Default is 5.
#' @param itermax Maximum number of iterations for k-means. Default is 100.
#' @param startn Number of random starts for k-means. Default is 10.
#'
#' @details
#' CASC integrates adjacency information with covariate similarity by forming
#' a combined matrix \eqn{Z Z^\top + \alpha C}, where \eqn{Z} is the normalized
#' adjacency and \eqn{C} is the covariate similarity. The tuning parameter
#' \eqn{\alpha} is selected by minimizing within-cluster sum of squares.
#'
#' @return A list containing:
#' \item{cluster}{Cluster membership assignments for each node.}
#' \item{alpha}{The selected tuning parameter value.}
#' \item{alphalower}{Lower bound of the search range for \eqn{\alpha}.}
#' \item{alphaupper}{Upper bound of the search range for \eqn{\alpha}.}
#'
#' @references
#' Binkiewicz, N., Vogelstein, J. T., & Rohe, K. (2017).
#' Covariate-assisted spectral clustering. *Biometrika*, 104(2), 361–377.
#'
#' @examples
#' # Example with toy adjacency and covariate matrices
#' set.seed(123)
#' Adj <- matrix(sample(0:1, 25, replace = TRUE), nrow = 5)
#' Adj[lower.tri(Adj)] <- t(Adj)[lower.tri(Adj)]  # make symmetric
#' diag(Adj) <- 0
#'
#' Covariate <- matrix(rnorm(5 * 2), nrow = 5, ncol = 2)
#'
#' casc_result <- CASC(Adj, Covariate, K = 2)
#' casc_result$cluster
#'
#' @export
CASC <- function (Adj, Covariate, K, alphan = 5, itermax = 100, startn = 10) {
  s = rowSums(Adj)
  s = s + mean(s)
  s = s^(-1/2)
  S = diag(s)
  Z = S %*% Adj %*% S
  net.eigen = eigen(Z %*% Z)
  ca = Covariate %*% t(Covariate)
  ca.eigen = eigen(ca)
  R <- ncol(Covariate)
  if (R <= K) {
    alphaupper <- net.eigen$values[1] / ca.eigen$values[R]
  } else {
    alphaupper <- net.eigen$values[1] / (ca.eigen$values[K] - ca.eigen$values[K + 1])
  }
  alphalower = (net.eigen$values[K] - net.eigen$values[K + 1]) / ca.eigen$values[1]
  d = rep(0, alphan)
  alpha = seq(alphalower, alphaupper, length.out = alphan)
  est = matrix(0, alphan, dim(Adj)[1])
  Norm <- function(x) sqrt(sum(x^2))
  for (ii in 1:alphan) {
    casc.eigen = eigen(Z %*% Z + alpha[ii] * ca)
    U = casc.eigen$vectors[, 1:K]
    Unorm = apply(U, 1, Norm)
    indu = which(Unorm > 0)
    U = U[indu, ] / Unorm[indu]
    result = kmeans(U, K, iter.max = itermax, nstart = startn)
    d[ii] = result$tot.withinss
    est[ii, indu] = as.factor(result$cluster)
  }
  best_index = which.min(d)
  result = est[best_index, ]
  best_alpha = alpha[best_index]
  return(list(cluster = result, alpha = best_alpha, alphalower = alphalower, alphaupper = alphaupper))
}
