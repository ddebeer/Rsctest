# Functions that are used in the process of computing the score-based tests
#  - do_pre_analysis: does some checking and handling of the arguments
#  - get_index: creates index based on order_by
#  - get_theta: estimates theta-values (using the 'PP'-package)
#  - get_stats: get the statistics for all order_by (i.e., indexes)
#  - get_stat: get the statistics for all which_col (for one order_by)


# function to create indexes based on list of order_by
get_index <- function(order_by, nPerson){
  nOrder <- length(order_by)
  index <- matrix(NA, nrow = nPerson, ncol = length(order_by),
                  dimnames = list(NULL, names(order_by)))
  for(orderNr in seq_len(nOrder)){
    index[,orderNr] <- order(order_by[[orderNr]])
  }
  return(index)
}

# function that returns theta (and checks its dim)
# or, when theta == NULL, estimates the theta-values using the pp-package
get_theta <- function(resp,
                      a = rep(1, length(b)),
                      b,
                      c = rep(0, length(b)),
                      theta = NULL,
                      theta_method = c("wle", "mle", "eap", "map"),
                      slope_intercept = FALSE){

  if(is.null(theta)){
    # compute person parameter estimates for all persons
    theta_method <- match.arg(theta_method)
    d <- rep(1, length(b))
    # check with IRT model formulation:
    if(slope_intercept){
      print("Using the 'PP'-package to estimate the ability parameters.")
      thetaEst <- PP::PP_4pl(resp, thres = b/a, slopes = a, lowerA = c, upperA = d,
                         type = theta_method)
    } else {
      thetaEst <- PP::PP_4pl(resp, thres = b, slopes = a, lowerA = c, upperA = d,
                         type = theta_method)
    }

    theta <- thetaEst$resPP$resPP[,1]
  } else if(length(theta) != dim(resp)[1])
    stop("'theta' should be a vector of length equal to the number of rows in 'resp'")

  return(theta)
}


# function to get test statistics for each order_by and each set of which_col
# returns a matrix with dimensions = c(# order_by, # which_col)
get_stats <- function(process, index, which_col, type){
  nIndex <- dim(index)[2]
  nWhich_col <- length(which_col)

  # create matrix wiht stats to return
  out <- matrix(NA, nrow = nIndex, ncol = nWhich_col,
                dimnames = list(colnames(index), names(which_col)))
  for(indexNr in seq_len(nIndex)){
    out[indexNr, ] <- get_stat(process, index[, indexNr],
                               type, which_col)
  }
  return(out)
}




# Function to compute the test statistic
get_stat <- function(process, index,
                     type = c("max", "doubleMax", "meanL2", "maxL2"),
                     which_col){

  # Order score contributions
  process <- process[index, , drop = FALSE]

  # compute cumSums
  B <- apply(process, 2, cumsum)

  # Compute Statistic
  nPer <- nrow(process)
  i_n <- (1:nPer) / nPer
  stat <- switch(type,
                 "max" = apply(B[,which_col], 2, function(x) max(abs(x))),
                 "doubleMax" = sapply(which_col, function(colNrs) {
                   max(apply(B[,colNrs], 2, function(x) max(abs(x))))
                 }),
                 "meanL2" = sapply(which_col, function(colNrs) {
                   mean(rowSums(B^2))
                 }),
                 "maxL2" = sapply(which_col, function(colNrs) {
                   max(( (i_n * (1 - i_n))^(-1) * rowSums(B^2)  ))
                 }))

  return(stat)
}

