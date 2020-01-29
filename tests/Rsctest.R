set.seed(7485)

# function to simulate 3PL data
P3pl <- function(theta, a = rep(1, length(b)), b, c = rep(0, length(b)),
                 generateResponses = FALSE, seed = 1){
  nItem <- length(b)
  nPerson <- length(theta)
  A <- rep(1, nPerson) %o% a
  B <- rep(1, nPerson) %o% b
  C <- rep(1, nPerson) %o% c
  Theta <- theta %o% rep(1, nItem)

  Theta_B <- Theta - B
  P <- C + (1 - C) * (1 + exp(- A * Theta_B))^(-1)
  if(generateResponses){
    set.seed(seed)
    resp <- (P > matrix(runif(nPerson * nItem),
                        ncol = nItem, nrow = nPerson)) * 1
    return(resp)
  } else (return(P))
}

# item parameters
b <- seq(-2, 1.8, length.out = 3)
a <- runif(length(b), 1, 1.5)
c <- rep(0, length(b))

# person parameters
nPerson <- 20
theta <- rnorm(nPerson)

# generate responses
resp <- P3pl(theta, a, b, generateResponses = TRUE)

# create orders
metric <- rnorm(nPerson)
factor <- factor(sample(1:3, size = nPerson, replace = TRUE))
order <- ordered(sample(1:3, size = nPerson, replace = TRUE))

# compute scores
library(Rsctest)
scores <- Rsctest:::get_scores(resp, a, b, c, theta,
                     slope_intercept = FALSE, return_terms = TRUE)
scores

# compute process
process <- apply(Rsctest:::scale_scores(scores$scores,
                                        decorrelate = TRUE)[order(metric),],
                 2, cumsum)
process

# Compare with strucchange
library(strucchange)
gefp <- gefp(x = scale(scores$scores, scale = FALSE), fit = NULL,
             scores = function(x) {x}, order.by = metric)

# Is the cumsum score process the same?
stopifnot(all((gefp$process[-1,] - process) < 1e-14))

# bootstrap vs sctest
test_b <- bootstrap_sctest(resp = resp, a = a, b = b, nSample = 15,
                         item_selection = 1,
                         DIF_covariate = list(m1 = metric,
                                         m2 = metric,
                                         m3 = metric,
                                         f = factor,
                                         o1 = order,
                                         o2 = order),
                         statistic = c("auto", "CvM", "maxLM", "auto", "auto", "WDMo"),
                         decorrelate = TRUE, theta = theta)
test_b$p
test_b$statistic
sctest_m1 <- sctest(x = scale(scores$scores, scale = FALSE),
                 scores = function(x) {x}, parm = c(1, 4), order.by = metric)
sctest_m2 <- sctest(x = scale(scores$scores, scale = FALSE),
                    scores = function(x) {x}, parm = c(1, 4), order.by = metric,
                    functional = "CvM")
sctest_m3 <- sctest(x = scale(scores$scores, scale = FALSE),
                    scores = function(x) {x}, parm = c(1, 4), order.by = metric,
                    functional = "maxLM")
sctest_f <- sctest(x = scale(scores$scores, scale = FALSE),
                   scores = function(x) {x}, parm = c(1, 4), order.by = factor,
                   functional = "LMuo")
sctest_o <- sctest(x = scale(scores$scores, scale = FALSE),
                   scores = function(x) {x}, parm = c(1, 4), order.by = order,
                   functional = "maxLMo")
sctest_o2 <- sctest(x = scale(scores$scores, scale = FALSE),
                   scores = function(x) {x}, parm = c(1, 4), order.by = order,
                   functional = "WDMo")

# Are the test statistics the same
stopifnot(all(test_b$statistic / rbind(sctest_m1$statistic,
                                      sctest_m2$statistic,
                                      sctest_m3$statistic,
                                      sctest_f$statistic,
                                      sctest_o$statistic,
                                      sctest_o2$statistic) - 1 < 1e-14))

# permutation vs bootstrap
test_p <- permutation_sctest(resp = resp, a = a, b = b, nSample = 15,
                             item_selection = 1,
                             DIF_covariate = list(m1 = metric,
                                             m2 = metric,
                                             m3 = metric,
                                             f = factor,
                                             o1 = order,
                                             o2 = order),
                             statistic = c("auto", "CvM", "maxLM", "auto", "auto", "WDMo"),
                             decorrelate = TRUE, theta = theta)
test_p$p
test_p$statistic

# Are the test statistics the same
stopifnot(all(test_b$statistic == test_p$statistic))


