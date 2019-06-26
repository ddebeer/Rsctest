#' A score-based DIF test using the permutation approach.
#'
#' \code{permutation_sctest} computes score test to detect DIF in multiple
#' item/parameters with respect to multiple person covariates (\code{order_by}).
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
#'   \item{\code{test_stats}}{A matrix containing all the test statistics.}
#'   \item{\code{p}}{A matrix containing the obtained \emph{p}-values.}
#'   \item{\code{nSamples}}{The number of samples taken.}
#'   \item{\code{order_by}}{A list containing all the covariate(s) used to order
#'    the score contirbutions.}
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
                               order_by = NULL,
                               parameters = c("per_item", "ab", "a", "b"),
                               itemNrs = NULL,
                               nSamples = 1000,
                               theta_method = c("wle", "mle", "eap", "map"),
                               slope_intercept = FALSE,
                               type = c("auto", "doubleMax", "meanL2", "maxL2"),
                               meanCenter = TRUE,
                               decorrelate = TRUE,
                               impact_groups = rep(1, dim(resp)[1])){


  # The response should be in a matrix
  stopifnot(is.matrix(resp) | is.data.frame(resp))
  if(is.data.frame(resp)) resp <- as.matrix(resp)

  # number of items, number of persons
  nItem <- ncol(resp)
  nPerson <- nrow(resp)

  # length of item parameter vectors should be correct
  if(!all(sapply(list(a, b, c), length) == nItem))
    stop("The vector of item parameters should be equal to the number of items in the respons matrix")

  # retrieve theta (or estimate when theta == NULL)
  theta <- get_theta(resp, a, b, c, theta, theta_method, slope_intercept)

  # match arguments
  type <- match.arg(type)
  parameters <- match.arg(parameters)

  # get list of which_col
  # select the items
  if(is.null(itemNrs)) itemNrs <- seq_len(nItem)
  if(!all(itemNrs %in% seq_len(nItem))) stop("'itemNrs' does not correspond with the number of items in the response matrix.")

  # select the parameters
  if(parameters == "per_item"){
    which_col <- lapply(itemNrs, function(itemNr) {
      c(itemNr, itemNr + nItem)
    })
    names(which_col) <- paste0("item", itemNrs)
    if(type == "auto") type <- "doubleMax"
  } else {
    which_col <- switch(parameters,
                        "a" = itemNrs,
                        "b" = itemNrs + nItem,
                        "ab" = c(itemNrs, itemNrs + nItem))
    names(which_col) <- switch(parameters,
                               "a" = paste0("a_it", itemNrs),
                               "b" = paste0("b_it", itemNrs),
                               "ab" = c(paste0("a_it", itemNrs),
                                        paste0("b_it", itemNrs)))

    if(type == "auto") type <- "max"
  }

  # If no order variable is given create one
  if(is.null(order_by)) order_by <- list(order1 = seq_len(nPerson))
  if(length(order_by) == nPerson) order_by <- list(order1 = order_by)

  # order_by should be a list with elements of length nPerson
  stopifnot(is.list(order_by))
  if(!all(sapply(order_by, function(x) length(x) == nPerson)))
    stop("The length of the element(s) of order_by is not equal to the number
of score contributions")

  # give names to order_by
  if(is.null(names(order_by)))
    names(order_by) <- paste0("order", seq_along(order_by))

  # create index- matrix according to the order_bys
  index <- get_index(order_by, nPerson)

  # get the scores, as well as the terms to compute the scores
  scores_terms <- get_scores(resp, a, b, c, theta,
                             slope_intercept, sparse = FALSE,
                             return_terms = TRUE)

  # scale generated score contributions, so that they become brownian process
  process <- scale_scores(scores_terms$scores, meanCenter, decorrelate, impact_groups)

  # compute the test statistic based on the observed scores
  test_stats <- get_stats(process, index, which_col, type)

  # get test statistic distribution based on the permutations
  permuted_stats <- get_permuted_stats(process, which_col, type, nSamples)

  # compute the p-values
  p <- get_pvalues(test_stats, permuted_stats)


  return(list(test_stats = test_stats,
              p = p,
              nSamples = nSamples,
              order_by = order_by))


}


# function to compute the bootstrapped statistics
get_permuted_stats <- function(process, which_col, type, nSamples){

  sapply(seq_len(nSamples), get_one_permuted_stat, simplify = "array",
         process, which_col, type)
}


# function to compute one bootstrapped statistic
get_one_permuted_stat <- function(sampleNr, process, which_col, type){

  nPerson <- 'if'(is.null(dim(process)), length(process), dim(process)[1])

  # permute the index
  index <- matrix(sample.int(nPerson, replace = FALSE), ncol =1)

  # compute statistics
  stats <- get_stats(process, index, which_col, type)

  return(stats)
}
