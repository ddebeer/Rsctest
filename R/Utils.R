# Utility functions that are used in the process of computing the score-based tests
#  - get_which_col: gets the columns for the tests
#  - get_theta: estimates theta-values (using the 'PP'-package)

# Function to compute the list with columnnumbers for the score tests
get_which_col <- function(item_selection, resp, parameters){
  
  # get number of items
  nItem <- dim(resp)[2]

  # Collect the column names
  colnames(resp) <- itemNames <- 'if'(is.null(colnames(resp)),
                                      sprintf(paste("it%0", nchar(nItem),
                                                    "d", sep=''),
                                              seq_len(nItem)),
                                      colnames(resp))



  # select the items
  itemNrs <-
    if(is.null(item_selection)) {seq_len(nItem)
    } else if(is.character(item_selection)){
      if(!all(item_selection %in% itemNames)) stop(
        "Some of the strings in 'item_selection'",
        " do not correspond to the column names of 'resp'.")
      else which(itemNames %in% item_selection)
    } else if(!all(item_selection %in% seq_len(nItem))) {stop(
      "Some values in'item_selection'",
      " do not correspond to the number of items in 'resp'.")
    } else item_selection

  if(parameters == "per_item"){

    which_col <- lapply(itemNrs, function(itemNr) {
      c(itemNr, itemNr + nItem)
    })
    names(which_col) <- itemNames[itemNrs]
  } else {
    which_col <- switch(parameters,
                        "a" = itemNrs,
                        "b" = itemNrs + nItem,
                        "ab" = c(itemNrs, itemNrs + nItem))
    names(which_col) <- switch(parameters,
                               "a" = paste0("a_", itemNames[itemNrs]),
                               "b" = paste0("b_", itemNames[itemNrs]),
                               "ab" = c(paste0("a_", itemNames[itemNrs]),
                                        paste0("b_", itemNames[itemNrs])))

  }
  return(which_col)
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
