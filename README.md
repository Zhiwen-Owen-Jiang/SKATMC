# SKATMC
The sequence kernel association test for multi-categorical outcomes (SKATMC)

## Installation
```{r}
# install.packages("devtools") if you have not installed devtools
devtools::install_github("Zhiwen-Owen-Jiang/SKATMC")
```
## Usage
```{r}
library(SKATMC)

# let's first generate some data
set.seed(123)
test.data <- data.frame(subtype = as.factor(sample(4, 100, replace = T)),
                              age = rnorm(n = 100, mean = 30, sd = 5),
                              sex = rbinom(100, 1, 1/2))
G <- matrix(rbinom(1000, 2, 0.05), nrow = 100, ncol = 10)

# Then do the analysis
test.null.model.glm <- SKATMC_Null_Model(subtype ~ age + sex, data.type = 'nominal',
                                     data = test.data, ref = 4) # generate the null model and set the 4th category as the reference
SKATMC(test.null.model.glm, G = G, weights = F) # do not use weight as in common variant analysis
SKATMC(test.null.model.glm, G = G, weights = T) # use the default weight as in rare variant analysis
SKATMC(test.null.model.glm, G = G, weights = T, weights.beta = c(1, 25)) # use specified weight, where c(1, 25) are shape parameters for the density of beta distribution
```
## Getting help
Please email Owen Jiang <owenjf@live.unc.edu>.

## Citation
Jiang, Zhiwen, et al. "The sequence kernel association test for multicategorical outcomes." Genetic epidemiology 47.6 (2023): 432-449.
