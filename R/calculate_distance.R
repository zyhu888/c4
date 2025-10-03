#' Calculate Distance Matrices from Covariates
#'
#' Compute pairwise distance matrices for a given dataset using one of three methods:
#' Euclidean, Hamming, or Gower distance. The function supports both numeric and
#' categorical variables, and automatically converts character columns to factors.
#'
#' @importFrom stats dist as.dist
#' @importFrom cluster daisy
#'
#' @param data A data frame or matrix containing covariates.
#' Character variables will be converted to factors automatically.
#' @param method A character string specifying the distance metric.
#' One of \code{"euclidean"}, \code{"hamming"}, or \code{"gower"}.
#'
#' @details
#' - \strong{Euclidean:} Requires all columns to be numeric.
#' - \strong{Hamming:} Accepts binary, logical, factor, or numeric-coded categorical data.
#' - \strong{Gower:} Requires the \pkg{cluster} package; handles mixed data types.
#'
#' @return A numeric distance matrix (class \code{matrix}).
#'
#' @examples
#' # Numeric example (Euclidean)
#' df_num <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' calculate_distance(df_num, method = "euclidean")
#'
#' # Categorical example (Hamming)
#' df_cat <- data.frame(a = c("A", "B", "A"), b = c("X", "X", "Y"))
#' calculate_distance(df_cat, method = "hamming")
#'
#' # Mixed data example (Gower)
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   df_mix <- data.frame(age = c(20, 25, 30),
#'                        gender = factor(c("M", "F", "M")))
#'   calculate_distance(df_mix, method = "gower")
#' }
#'
#' @export
calculate_distance <- function(data, method = c("euclidean", "hamming", "gower")) {
  method <- match.arg(method)

  if (!is.data.frame(data)) {
    stop("Input must be a data frame or matrix.")
  }

  # Convert character to factor
  data[] <- lapply(data, function(x) if (is.character(x)) as.factor(x) else x)

  if (method == "euclidean") {
    # Ensure all columns are numeric
    if (!all(sapply(data, is.numeric))) {
      stop("Euclidean distance requires all numeric columns.")
    }
    dist_matrix <- dist(data, method = "euclidean")

  } else if (method == "hamming") {
    # Convert to matrix for logical comparison
    if (!all(sapply(data, function(x) is.factor(x) || is.logical(x) || is.numeric(x)))) {
      stop("Hamming distance requires binary or categorical data.")
    }

    data_mat <- as.matrix(data)
    n <- nrow(data_mat)
    dist_matrix <- matrix(0, n, n)

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        dist_matrix[i, j] <- sum(data_mat[i, ] != data_mat[j, ])
        dist_matrix[j, i] <- dist_matrix[i, j]
      }
    }

    dist_matrix <- as.dist(dist_matrix)

  } else if (method == "gower") {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Please install the 'cluster' package to use Gower distance.")
    }
    dist_matrix <- cluster::daisy(data, metric = "gower")
  }

  return(as.matrix(dist_matrix))
}
