#' Generate a Weighted Stochastic Block Model (SBM) Network
#'
#' Simulate an undirected weighted network from a stochastic block model (SBM).
#' Within-community and between-community edges are drawn with specified
#' probabilities and weight distributions.
#'
#'@importFrom stats rbinom
#'
#' @param node_num Total number of nodes. Default is 200.
#' @param cluster_size A vector of community sizes. Must sum to \code{node_num}.
#' @param win_cluster_den Probability of edge existence within a community. Default is 0.3.
#' @param win_cluster_dist Distribution function for within-community edge weights (as string). Default is \code{"runif"}.
#' @param win_cluster_par List of parameters for \code{win_cluster_dist}. Default is \code{list(min = 1, max = 2)}.
#' @param btw_cluster_den Probability of edge existence between communities. Default is 0.1.
#' @param btw_cluster_dist Distribution function for between-community edge weights (as string). Default is \code{"runif"}.
#' @param btw_cluster_par List of parameters for \code{btw_cluster_dist}. Default is \code{list(min = 1, max = 2)}.
#'
#' @return A list containing:
#' \item{node_number}{Total number of nodes.}
#' \item{cluster_number}{Number of clusters.}
#' \item{adjacency_matrix}{Symmetric weighted adjacency matrix.}
#' \item{network_sparsity}{Proportion of absent edges in the network.}
#'
#' @examples
#' set.seed(123)
#' net <- gen_weighted_sbm(
#'   node_num = 20,
#'   cluster_size = c(10, 10),
#'   win_cluster_den = 0.5,
#'   btw_cluster_den = 0.1
#' )
#'  W <- net$adjacency_matrix
#'  W
#'
#' @export
gen_weighted_sbm <- function(node_num = 200,
                             cluster_size = c(50, 50, 50, 50),
                             win_cluster_den = 0.3,
                             win_cluster_dist = "runif",
                             win_cluster_par = list(min = 1, max = 2),
                             btw_cluster_den = 0.1,
                             btw_cluster_dist = "runif",
                             btw_cluster_par = list(min = 1, max = 2)) {
  if (node_num != sum(cluster_size)) {
    stop("The number of nodes is not equal to the sum of cluster sizes!")
  }

  cluster_number <- length(cluster_size)
  if (cluster_number < 2) {
    warning("The number of clusters is less than 2!")
  }

  node_member <- unlist(lapply(1:cluster_number, function(x) {
    rep(x, cluster_size[x])
  }))

  adj <- matrix(NA, nrow = node_num, ncol = node_num)

  for (i in 1:(node_num - 1)) {
    adj[i, i] <- 0
    for (j in (i + 1):node_num) {
      if (node_member[i] == node_member[j]) {
        adj[i, j] <- rbinom(1, 1, win_cluster_den) *
          do.call(win_cluster_dist, c(list(n = 1), win_cluster_par))
      } else {
        adj[i, j] <- rbinom(1, 1, btw_cluster_den) *
          do.call(btw_cluster_dist, c(list(n = 1), btw_cluster_par))
      }
      adj[j, i] <- adj[i, j]
    }
  }
  adj[node_num, node_num] <- 0

  if (!all(adj == t(adj))) {
    adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]
  }

  net_sparse <- 1 - sum(adj != 0) / (node_num * (node_num - 1) / 2)

  return(list(
    node_number = node_num,
    cluster_number = cluster_number,
    adjacency_matrix = adj,
    network_sparsity = net_sparse
  ))
}
