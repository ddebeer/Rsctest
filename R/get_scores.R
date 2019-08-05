# function to compute scores based on
# - resp = the response matrix (or data.frame)
# - a = a vector with the 'true' slope parameters (3PLM)
# - b = a vector with the 'true' location parameters (3PLM)
# - c = a vector with the 'true' pseudo-guessing parameters (3PLM)
# - d = a vector with the 'true' upper-limit parameters (3PLM)
# - theta = a vector with the 'true' ability parameters
# - slope_intercept = a logical specifying if the slope-intercept
#                    parameterization should be used.
#                    (FALSE = IRT-parameterization)
# - sparse = a logical indicating if a sparse matrix of the score-contributions
#            should be returned.
# - return_terms = a logical indicationg if the terms should be returned

get_scores <- function(resp, a, b, c, theta,
                       slope_intercept = FALSE,
                       sparse = FALSE,
                       return_terms = FALSE){

  # number of items, number of persons
  nItem <- ncol(resp)
  nPerson <- nrow(resp)

  # compute terms for scores
  terms <- get_terms_for_scores(theta, a, b, c, slope_intercept)

  # compute scores based on responses and terms
  scores <- get_scores_from_terms(resp, terms)

  if(sparse) scores <- Matrix::Matrix(scores)

  if(return_terms){
    return(list(scores = scores,
                terms = terms))
  }

  return(scores)
}


# function to compute the terms needed for computing the score contributions
get_terms_for_scores <- function(theta, a, b, c, slope_intercept){
  # number of items, number of persons
  nItem <- length(b)
  nPerson <- length(theta)

  # set item parameters and person parameter in matrix form
  A <- rep(1, nPerson) %o% a
  B <- rep(1, nPerson) %o% b
  C <- rep(1, nPerson) %o% c
  Theta <- theta %o% rep(1, nItem)

  # compute response probabilities
  if(slope_intercept){
    P <- C + (1 - C) * (1 + exp(- A * Theta + B))^(-1)
  } else {
    Theta_B <- Theta - B
    P <- C + (1 - C) * (1 + exp(- A * Theta_B))^(-1)
  }

  # compute P * (1 - C)
  P1_C <- P * (1 - C)

  # compute (P - C) / P * (1 - C)
  PC_P1C <- (P - C) / P1_C

  # compute terms to allow for easy score computation
  if(slope_intercept){
    term_a <- Theta * PC_P1C      # (Baker and Kim, 2004, p.47)
    term_b <- - PC_P1C            # (Baker and Kim, 2004, p.48)
  } else {
    term_a <- Theta_B * PC_P1C    # (Baker and Kim, 2004, p.47)
    term_b <- - A * PC_P1C        # (Baker and Kim, 2004, p.48)
  }
  term_c <- 1 / P1_C              # (Baker and Kim, 2004, p.48)

  return(list(P = P,
              a = term_a,
              b = term_b,
              c = term_c))
}


# function to compute the score contributions based on the terms
get_scores_from_terms <- function(resp, terms){
  X_P <- resp - terms$P

  score_a <- X_P * terms$a    # (Baker and Kim, 2004, p.47)
  score_b <- X_P * terms$b    # (Baker and Kim, 2004, p.48)
  score_c <- X_P * terms$c    # (Baker and Kim, 2004, p.48)

  # Calculating and returning the scores for a and b
  scores <- cbind(score_a,
                  score_b
                  # , score_c
                  )

  return(scores)

}


# function to scale/decorelate the scores (per impactGroup)
scale_scores <- function(scores, meanCenter = TRUE, decorrelate = TRUE,
                         impact_groups = rep(1, dim(scores)[1])){

  # create return object
  process <- scores

  for(groupIndicator in unique(impact_groups)){
    which <- impact_groups == groupIndicator

    if(decorrelate){
      # Center scores
      process[which, ] <- 'if'(meanCenter, scale(process[which, ], scale = FALSE), process[which, ])

      # missing response (NA) corresponds with zero (0) score contribution
      process[which, ][is.na(process[which, ])] <- 0

      # Decorelate process
      sqrtInvI <- tryCatch(
        chol2inv(chol(expm::sqrtm(crossprod(process[which, ])))),
        error = stop("Decorrelating leads to numerical problems.
Use 'decorrelate = FALSE' instead.",
                     call. = FALSE)
        )
      process[which, ] <- (process[which, ] %*% sqrtInvI)
    } else {
      # Center scores
      process[which, ] <- scale(process[which, ], center = meanCenter)
      # missing response (NA) corresponds with zero (0) score contribution
      process[which, ][is.na(process[which, ])] <- 0
    }
  }

  return(process)
}

