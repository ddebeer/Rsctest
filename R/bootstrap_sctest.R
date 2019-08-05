#' A score-based DIF test using the parameteric bootstrap approach.
#'
#' \code{bootstrap_sctest} computes score test to detect DIF in multiple
#' item/parameters with respect to multiple person covariates (\code{order_by}).
#' To obtain the p-values a resampling approach is applied. That is, given the
#' (item and person) parameters, new data sets are sampled to create the
#' distribution of the test statistic under the null-hypothesis. The
#' functionality is limited to the 1-, 2-, and 3-parameter logistic models.
#' Only DIF with respect to the \code{a} and \code{b} parameters are tested for,
#' respectively the item discriminations and the item difficulties.
#'
#' @param resp A matrix (or data frame) containing the responses, with the
#' items in the columns.
#' @param a A vector of item slopes/item discriminations.
#' @param b A vector of item locations/item difficulties.
#' @param c A vector of pseudo guessing parameters.
#' @param theta A vector with the true/estimated ability parameters or NULL
#' (the default) which leads to the ability parameters being estimated.
#' @param order_by A list with the person covariate(s) to test for as
#' element(s).
#' @param parameters A charachter string, either "per_item", "ab", "a", or "b",
#' to specify which parameters should be tested for.
#' @param itemNrs An integer vector with the colum numbers in the \code{resp},
#' specifying the items for which the test should be computed. Or NULL (the
#' default), which leads to all the items being tested.
#' @param nSamples An integer value with the number of permutations to be
#' sampled.
#' @param theta_method A charachter string, either "wle", "mle", "eap", of
#' "map" that specifies the estimator for the ability estimation. Only
#' relevant when \code{theta == NULL}.
#' @param slope_intercept A logical value indicating whether the slope-intercept
#' formulation of the 2-/3-PL model should be used.
#' @param type A charachter string, either "auto", "doubleMax", "meanL2",
#' or "maxL2", specifying the test statistic to be used.
#' @param meanCenter A logical value: should the score contributions be mean
#' centered per parameter?
#' @param decorrelate A logical value: should the score contributions be
#' decorrelated?
#' @param impact_groups A vector indicating impact-group membership for
#' each person.
#' @return a list with four elements:
#' \describe{
#'   \item{\code{test_stats}}{A matrix containing all the test statistics.}
#'   \item{\code{p}}{A matrix containing the obtained \emph{p}-values.}
#'   \item{\code{nSamples}}{The number of samples taken.}
#'   \item{\code{order_by}}{A list containing all the covariate(s) used to order
#'    the score contirbutions.}
#' }
#' @aliases bootstrap_sctest
#' @seealso \code{\link{permutation_sctest}}
#'
#' @export
bootstrap_sctest <- function(resp,
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
                             decorrelate = FALSE,
                             impact_groups = rep(1, dim(resp)[1])){


  # The responses should be in a matrix
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

  # get test statistic distribution based on
  bootstrapped_stats <- get_bootstrapped_stats(observed_resp = resp,
                                             terms = scores_terms$terms,
                                             meanCenter, decorrelate,
                                             impact_groups,
                                             index, which_col, type,
                                             nSamples)

  # compute the p-values
  p <- get_pvalues(test_stats, bootstrapped_stats)


  return(list(test_stats = test_stats,
              p = p,
              nSamples = nSamples,
              order_by = order_by,
              theta = theta))


}


# function to compute the bootstrapped statistics
get_bootstrapped_stats <- function(observed_resp, terms, meanCenter, decorrelate,
                                  impact_groups, index, which_col, type,
                                  nSamples = 1000
                                  ){

  # cat("Generating bootstrapped samples:\n")
  # progressBar <- txtProgressBar(min = 0, max = nSamples, style = 3, char = "|")
  # env <- environment()

  bootstrapped_stats <- lapply(
    seq_len(nSamples), get_one_bootstrapped_stat, observed_resp, terms,
    meanCenter, decorrelate, impact_groups, index, which_col, type
    #, env
  )

  array(unlist(bootstrapped_stats),
        dim = c(dim(bootstrapped_stats[[1]]), nSamples))
}


# function to compute one bootstrapped statistic
get_one_bootstrapped_stat <- function(sampleNr, observed_resp, terms, meanCenter, decorrelate,
                                      impact_groups, index, which_col, type
                                      #, env
                                      ){

  # setTxtProgressBar(get("progressBar", env), sampleNr)
  # generate response matrix
  gen_resp <- generate_response_matrix(terms$P, observed_resp)

  # compute scores based no generated response matrix
  gen_scores <- get_scores_from_terms(gen_resp, terms)

  # scale generated score contributions, so that they become brownian process
  process <- scale_scores(gen_scores, meanCenter, decorrelate, impact_groups)

  # compute statistics
  stats <- get_stats(process, index, which_col, type)

  return(stats)
}


# function to generate response matrix based on response probabilities (P)
# with the same pattern of missing responses as the oberved responses
generate_response_matrix <- function(P, observed_resp){
  # generate responses
  dims <- dim(observed_resp)
  generated_resp <- (P > stats::runif(prod(dims))) * 1
  dim(generated_resp) <- dims

  # missing in real data should be missing in generated data
  generated_resp[is.na(observed_resp)] <- NA

  return(generated_resp)
}


# function to get the p-values based on the computed test statistics and
# the bootstrapped/permuted statistics
get_pvalues <- function(test_stats, sampled_stats){
  p <- array(NA, dim = dim(test_stats))

  for(orderNr in seq_len(dim(p)[1])){
    for(which_colNr in seq_len(dim(p)[2])){
      p[orderNr, which_colNr] <- mean(sampled_stats[orderNr, which_colNr, ] >=
                                        test_stats[orderNr, which_colNr])
    }
  }

  dimnames(p) <- dimnames(test_stats)
  return(p)

}

