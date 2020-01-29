#' A score-based DIF test using the permutation approach.
#'
#' \code{permutation_sctest} computes score test to detect DIF in multiple
#' item/parameters with respect to multiple person covariates (\code{DIF_covariate}).
#' To obtain the p-values a resampling approach is applied. That is, person
#' orders are randomly permuted to sample from the test statistic distribution
#' under the null hypothesis. The
#' functionality is limited to the 1-, 2-, and 3-parameter logistic models.
#' Only DIF with respect to the \code{a} and \code{b} parameters are tested for,
#' respectively the item discriminations and the item difficulties.
#'
#' @inheritParams bootstrap_sctest
#'
#' @return a list with four elements:
#' \describe{
#'   \item{\code{statistics}}{A matrix containing all the test statistics.}
#'   \item{\code{p}}{A matrix containing the obtained \emph{p}-values.}
#'   \item{\code{nSamples}}{The number of samples taken.}
#'   \item{\code{DIF_covariate}}{A list containing all the covariate(s) used to order
#'    the score contirbutions, as well as the used test statistics.}
#' }
#' @aliases permutation_sctest
#' @seealso \code{\link{bootstrap_sctest}}
#'
#' @export
permutation_sctest <- function(resp,
                               theta = NULL,
                               a = rep(1, length(b)),
                               b,
                               c = rep(0, length(b)),
                               DIF_covariate = NULL,
                               parameters = c("per_item", "ab", "a", "b"),
                               item_selection = NULL,
                               nSamples = 1000,
                               theta_method = c("wle", "mle", "eap", "map"),
                               slope_intercept = FALSE,
                               statistic = "auto",
                               meanCenter = TRUE,
                               decorrelate = FALSE,
                               impact_groups = rep(1, dim(resp)[1])){


  # get call
  call <- match.call()
  
  # The responses should be in a matrix
  stopifnot(is.matrix(resp) | is.data.frame(resp))
  if(is.data.frame(resp)) resp <- as.matrix(resp)

  # number of persons
  nPerson <- nrow(resp)

  # retrieve theta (or estimate when theta == NULL)
  theta <- get_theta(resp, a, b, c, theta, theta_method, slope_intercept)

  # get list of which_col
  which_col <- get_which_col(item_selection, resp,
                             parameters = match.arg(parameters))

  # create index- matrix according to the DIF_covariates
  index_list <- get_index_list(DIF_covariate, nPerson, statistic, call)

  # get the scores, as well as the terms to compute the scores
  scores_terms <- get_scores(resp, a, b, c, theta,
                             slope_intercept, sparse = FALSE,
                             return_terms = TRUE)

  # scale generated score contributions, rather than scaling the brownian process
  scaled_scores <- scale_scores(scores_terms$scores, meanCenter, decorrelate,
                                impact_groups)

  # compute the test statistic based on the observed scores
  test_stats <- get_stats(scaled_scores, index_list, which_col)

  # get test statistic distribution based on the permutations
  permuted_stats <- get_permuted_stats(scaled_scores, which_col,
                                       index_list, nSamples)

  # compute the p-values
  p <- get_pvalues(test_stats, permuted_stats)


  return(list(resp = resp, 
              statistic = test_stats,
              p = p,
              nSamples = nSamples,
              DIF_covariate = index_list,
              theta = theta))


}


# function to compute the bootstrapped statistics
get_permuted_stats <- function(scaled_scores, which_col, index_list, nSamples){

  permuted_stats <- lapply(
    seq_len(nSamples), get_one_permuted_stat, scaled_scores, which_col, index_list)

  array(unlist(permuted_stats),
        dim = c(dim(permuted_stats[[1]]), nSamples))
}


# function to compute one bootstrapped statistic
get_one_permuted_stat <- function(sampleNr, scaled_scores, which_col, index_list){

  nPerson <- 'if'(is.null(dim(scaled_scores)), length(scaled_scores), dim(scaled_scores)[1])

  # permute the index
  permuted_index <- sample.int(nPerson, replace = FALSE)

  # compute statistics
  stats <- get_stats(scaled_scores, index_list, which_col,
                     permuted_index = permuted_index)

  return(stats)
}
