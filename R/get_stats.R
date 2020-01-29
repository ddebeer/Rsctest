# Functions that are used in the process of computing the score-based tests
#  - get_index_list: creates list of index + statistic based on DIF_covariate
#  - create_statistic: function that creates the function to compute the test
#       statistic
#  - get_stats: get the statistics for all DIF_covariate (i.e., indexes)
#  - get_stat: get the statistics for all which_col (for one DIF_covariate)


# function to create a list of indexes based on list of DIF_covariate
get_index_list <- function(DIF_covariate, nPerson, statistic, call){

  # If no order variable is given create one
  if(is.null(DIF_covariate)) DIF_covariate <- list(order1 = seq_len(nPerson))
  if(length(DIF_covariate) == nPerson && length(DIF_covariate[[1]]) == 1){
    DIF_covariate <- list(DIF_covariate)
    names(DIF_covariate) <- as.character(deparse(call$DIF_covariate)) 
  }


  # DIF_covariate should be a list with elements of length nPerson
  stopifnot(is.list(DIF_covariate))
  if(!all(sapply(DIF_covariate, function(x) length(x) == nPerson)))
    stop("The length of the element(s) of DIF_covariate is not equal to the number
         of score contributions")

  # give names to DIF_covariate
  if(is.null(names(DIF_covariate)) & length(DIF_covariate) > 1) {
    # get names from call to previous function.
    names(DIF_covariate) <- as.character(match.call(sys.function(-1),
                                               sys.call(-1))$DIF_covariate[-1])
  }

  # check if the statistics are OK
  stopifnot(all(statistic %in% c("auto", "DM", "CvM", "maxLM", "LMuo",
                                 "WDMo", "maxLMo")))

  # check length of statistic
  nDIF_covariate <- length(DIF_covariate)
  if(length(statistic) == 1 & nDIF_covariate > 1)
    statistic <- rep(statistic, nDIF_covariate)
  stopifnot(length(statistic) == nDIF_covariate)

  # create index list with each object corresponding with one DIF_covariate
  index_list <- DIF_covariate
  for(orderNr in seq_len(nDIF_covariate)){
    this_DIF_covariate <- DIF_covariate[[orderNr]]
    this_index <- order(this_DIF_covariate)
    this_statistic <- statistic[orderNr]
    index_list[[orderNr]] <- list(
      index = this_index,
      DIF_covariate = this_DIF_covariate[this_index],
      statistic = create_statistic(this_DIF_covariate, this_statistic)
    )
    }

  return(index_list)
}


# function to create the test statistic, based on an DIF_covariate variable
create_statistic <- function(DIF_covariate, statistic){

  # Select the defaults
  if(statistic == "auto"){

    # get the variable type for the DIF_covariate variable
    type_DIF_covariate <- get_variable_type(DIF_covariate)

    statistic <- switch(type_DIF_covariate,
                        "metr" = "DM",
                        "cat" = "LMuo",
                        "ordcat" = "maxLMo")
  }

  # select the test statistic
  statistic <- switch(statistic,
                      "DM" = DM(),
                      "CvM" = CvM(),
                      "maxLM" = maxLM(DIF_covariate),
                      "LMuo" = LMuo(DIF_covariate),
                      "WDMo" = WDMo(DIF_covariate),
                      "maxLMo" = maxLMo(DIF_covariate)# ,
                      # stop("stat should be on of the following character strings: ",
                      #     "'auto', 'DM', 'CvM', 'maxLM', 'LMuo', or 'WDMo', 'maxLMo'")
                      )

  return(statistic)
}


# function that computes the maximum absolute value of a vector (with weigths)
wmax_abs <- function(x, weights = 1) suppressWarnings(max(abs(x * weights)))

# functions to create the statistics
### Double Max
DM <- function()
  list(stat = function(process, which_col)
    sapply(which_col, function(colNrs) {
      max(apply(process[,colNrs, drop = FALSE], 2, wmax_abs))
      }),
    name = "Double Maximum Test")

### Cramer-von Mises
CvM <- function()
  list(stat = function(process, which_col)
    sapply(which_col, function(colNrs) {
      mean(rowSums(process[,colNrs, drop = FALSE]^2))
    }),
    name = "Cramer-von Mises Test")

### Maximum Lagrange Multiplier Test
maxLM <- function(DIF_covariate){
  nPer <- length(DIF_covariate)
  i_n <- (1:nPer) / nPer
  weights <- (i_n * (1 - i_n))[-nPer]

  list(stat = function(process, which_col){
    sapply(which_col, function(colNrs) {
      max(( rowSums(process[-nPer,colNrs, drop = FALSE]^2) / weights))
    })
  },
  name = "Maximum Lagrange Multiplier Test")
}

### Maximum Lagrange Multiplier Test for Unordered Groups
LMuo <- function(DIF_covariate){
  freq <- table(DIF_covariate)
  prop <- table(DIF_covariate) / length(DIF_covariate)
  i_cat <- cumsum(freq)
  catdiffL2 <- function(column)
    sum(diff(c(0, column[i_cat]))^2 / freq * length(column))

  list(stat = function(process, which_col){
    sapply(which_col, function(colNrs) {
      sum(apply(process[,colNrs, drop = FALSE], 2, catdiffL2))
    })
  },
  name = "Lagrange Multiplier Test for Unordered Groups")
}

### Maximum Lagrange Multiplier Test for Ordered Groups
maxLMo <- function(DIF_covariate){
  freq <- table(DIF_covariate)
  i_cat <- cumsum(freq)[-length(freq)]
  i_n <- i_cat / length(DIF_covariate)
  weights <- (i_n * (1 - i_n))

  list(stat = function(process, which_col){
    sapply(which_col, function(colNrs) {
      max(( rowSums(process[i_cat,colNrs, drop = FALSE]^2) / weights))
    })
  }, name = "Maximum Lagrange Multiplier Test for Ordered Groups")
}

### Weighted Double Maximum for Ordered Groups
WDMo <- function(DIF_covariate){
  freq <- table(DIF_covariate)
  i_cat <- cumsum(freq)[-length(freq)]
  i_n <- i_cat / length(DIF_covariate)
  weights <- 1 / sqrt(i_n * (1 - i_n))

  list(stat = function(process, which_col){
    sapply(which_col, function(colNrs) {
      max(apply(process[i_cat,colNrs, drop = FALSE], 2, wmax_abs, weights = weights))
    })
  }, name = "Weighted Double Maximum for Ordered Groups")
}


# Function to get the variable type of a variable (here an DIF_covariate variable)
get_variable_type <- function(variable){
  class <- class(variable)
  if("ordered" %in% class) {
    "ordcat"
  } else if (class %in% c("factor", "logical", "character")){
    "cat"
  } else if (class %in% c("integer", "numeric")){
    "metr"
  } else stop("Change the class of the variable")

}

# function to get test statistics for each DIF_covariate and each set of which_col
# returns a matrix with dimensions = c(# DIF_covariate, # which_col)
get_stats <- function(scaled_scores, index_list, which_col, permuted_index = NULL){
  nIndex <- length(index_list)
  nWhich_col <- length(which_col)

  # create matrix wiht stats to return
  out <- matrix(NA, nrow = nIndex, ncol = nWhich_col,
                dimnames = list(names(index_list), names(which_col)))
  for(indexNr in seq_len(nIndex)){
    index <- 'if'(is.null(permuted_index),
                  index_list[[indexNr]]$index,
                  permuted_index)
    out[indexNr, ] <- get_stat(scaled_scores, index,
                               index_list[[indexNr]]$statistic$stat,
                               which_col)
  }
  return(out)
}




# Function to compute the test statistic
get_stat <- function(scaled_scores, index, stat, which_col){

  # Order score contributions
  scaled_scores <- scaled_scores[index, , drop = FALSE]

  # compute cumSums
  process <- apply(scaled_scores, 2, cumsum)

  # Compute Statistic
  stat <- stat(process, which_col)
  return(stat)
}

