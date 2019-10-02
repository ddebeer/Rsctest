# Functions that are used in the process of computing the score-based tests
#  - get_index_list: creates list of index + statistic based on order_by
#  - create_statistic: function that creates the function to compute the test
#       statistic
#  - get_stats: get the statistics for all order_by (i.e., indexes)
#  - get_stat: get the statistics for all which_col (for one order_by)


# function to create a list of indexes based on list of order_by
get_index_list <- function(order_by, nPerson, statistic){

  # If no order variable is given create one
  if(is.null(order_by)) order_by <- list(order1 = seq_len(nPerson))
  if(length(order_by) == nPerson && length(order_by[[1]]) == 1){
    order_by <- list(order_by)
    names(order_by) <- as.character(match.call(sys.function(-1),
                                               sys.call(-1))$order_by)
  }


  # order_by should be a list with elements of length nPerson
  stopifnot(is.list(order_by))
  if(!all(sapply(order_by, function(x) length(x) == nPerson)))
    stop("The length of the element(s) of order_by is not equal to the number
         of score contributions")

  # give names to order_by
  if(is.null(names(order_by)) & length(order_by) > 1) {
    # get names from call to previous function.
    names(order_by) <- as.character(match.call(sys.function(-1),
                                               sys.call(-1))$order_by[-1])
  }

  # check if the statistics are OK
  stopifnot(all(statistic %in% c("auto", "DM", "CvM", "maxLM", "LMuo",
                                 "WDMo", "maxLMo")))

  # check length of statistic
  nOrder_by <- length(order_by)
  if(length(statistic) == 1 & nOrder_by > 1)
    statistic <- rep(statistic, nOrder_by)
  stopifnot(length(statistic) == nOrder_by)

  # create index list with each object corresponding with one order_by
  index_list <- order_by
  for(orderNr in seq_len(nOrder_by)){
    this_order_by <- order_by[[orderNr]]
    this_index <- order(this_order_by)
    this_statistic <- statistic[orderNr]
    index_list[[orderNr]] <- list(
      index = this_index,
      order_by = this_order_by[this_index],
      statistic = create_statistic(this_order_by, this_statistic)
    )
    }

  return(index_list)
}


# function to create the test statistic, based on an order_by variable
create_statistic <- function(order_by, statistic){

  # Select the defaults
  if(statistic == "auto"){

    # get the variable type for the order_by variable
    type_order_by <- get_variable_type(order_by)

    statistic <- switch(type_order_by,
                        "metr" = "DM",
                        "cat" = "LMuo",
                        "ordcat" = "maxLMo")
  }

  # select the test statistic
  statistic <- switch(statistic,
                      "DM" = DM(),
                      "CvM" = CvM(),
                      "maxLM" = maxLM(order_by),
                      "LMuo" = LMuo(order_by),
                      "WDMo" = WDMo(order_by),
                      "maxLMo" = maxLMo(order_by)# ,
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
    name = "Double Maximum")

### Cramer-von Mises
CvM <- function()
  list(stat = function(process, which_col)
    sapply(which_col, function(colNrs) {
      mean(rowSums(process[,colNrs, drop = FALSE]^2))
    }),
    name = "Cramer-von Mises")

### Maximum Lagrange Multiplier Test
maxLM <- function(order_by){
  nPer <- length(order_by)
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
LMuo <- function(order_by){
  freq <- table(order_by)
  prop <- table(order_by) / length(order_by)
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
maxLMo <- function(order_by){
  freq <- table(order_by)
  i_cat <- cumsum(freq)[-length(freq)]
  i_n <- i_cat / length(order_by)
  weights <- (i_n * (1 - i_n))

  list(stat = function(process, which_col){
    sapply(which_col, function(colNrs) {
      max(( rowSums(process[i_cat,colNrs, drop = FALSE]^2) / weights))
    })
  }, name = "Maximum Lagrange Multiplier Test for Ordered Groups")
}

### Weighted Double Maximum for Ordered Groups
WDMo <- function(order_by){
  freq <- table(order_by)
  i_cat <- cumsum(freq)[-length(freq)]
  i_n <- i_cat / length(order_by)
  weights <- 1 / sqrt(i_n * (1 - i_n))

  list(stat = function(process, which_col){
    sapply(which_col, function(colNrs) {
      max(apply(process[i_cat,colNrs, drop = FALSE], 2, wmax_abs, weights = weights))
    })
  }, name = "Weighted Double Maximum for Ordered Groups")
}


# Function to get the variable type of a variable (here an order_by variable)
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

# function to get test statistics for each order_by and each set of which_col
# returns a matrix with dimensions = c(# order_by, # which_col)
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

