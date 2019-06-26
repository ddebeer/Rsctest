### Resampling Score(-based) Test for DIF Detection.

This R-package includes two approaches to perform a score-test for DIF detection based on resampling techniques, and is focused on data from Multi-Stage Adaptive Tests. Currently the functionality is limited to the 1-, 2-, and 3-PL models, and only the DIF in the threshold (b) and slope (a) parameter are tested for. Rather than using the theoretical distribution to obtain the __p__-value for the score-based test statistic, this distribution is obtained using resampling techniques. 

- `permutation_sctest` samples permutations of person orders.

- `bootstrap_sctest` samples new data based on the estimated (item and person).


## Installation


The package can be installed using using the `devtools`-package:

```
install.packages("devtools")
devtools::install_github("ddebeer/Rsctest")
```



