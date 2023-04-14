#' Get paramters of the NULL model
#'
#' @param formula an object of class "formula": a symbolic description of the
#' NULL model to be fitted.
#' @param data.type an indicator of the outcome type, either 'nominal' or 'ordinal'.
#' @param ref an optional argument for specifying the reference group, only
#' valid for nominal data. It can be a scalar object or a vector containing levels
#' of the outcome. The default is `NULL` where every level of the outcome will
#' be tried as the reference.
#' @param data an optional data frame containing the variables in the model
#' (default=NULL).
#'
#' @return
#' This function returns an object that has model parameters of
#' the NULL model of no association between genetic variables and outcomes.
#' After obtaining it, please use `SKATMC` function to conduct the association test.
#' @export
#'
#' @examples
SKATMC_Null_Model <- function(formula, data.type, ref = NULL, data = NULL){
  # browser()
  if (data.type == 'nominal'){
    invisible(utils::capture.output(fit <- mclogit::mblogit(formula, data, na.action = 'na.fail')))
    X <- model.matrix(formula, fit$model)
    nSam <- nrow(X)
    y.matrix <- matrix(fit$y, nrow = nSam, byrow = T)
    J <- ncol(y.matrix)
    mu.hat <- matrix(fit$fitted.values, ncol = J, byrow = T)

    if (is.null(ref)){
      all.levels <- rownames(fit$D)
    } else {
      if (!all(ref %in% rownames(fit$D))) stop('ref should be a level of outcome.')
      all.levels <- ref
    }

    len.all.levels <- length(all.levels)
    W.list <- y.star.list <- rep(0, len.all.levels)
    names(W.list) <- names(y.star.list) <- all.levels

    for (ref.i in all.levels){
      ref.index <- which(rownames(fit$D) == ref.i)
      D <- c((mu.hat[, -ref.index] + mu.hat[, ref.index]) / (mu.hat[, -ref.index] * mu.hat[, ref.index]))
      # D <- 1/D
      y.star <- list(D * (c(y.matrix[, -ref.index]) - c(mu.hat[, -ref.index])))
      y.star.list[as.character(ref.i)] <- y.star

      WI <- get_WI_GLM(mu.hat[, -ref.index], D = D, J = J, nSam = nSam)
      W <- list(quickInverse(WI, dimen = J))
      W.list[as.character(ref.i)] <- W
    }

  } else if (data.type == 'ordinal'){
    fit <- MASS::polr(formula, data, method = 'logistic', na.action = 'na.fail')
    if (is.null(data)) data <- fit$model
    X <- model.matrix(formula, data)
    nSam <- nrow(X)
    y <- as.formula(paste0('~ 0 +', as.character(formula)[[2]]))
    y.matrix <- model.matrix(y, data)
    J <- ncol(y.matrix)

    pi.hat <- fit$fitted.values
    A <- matrix(1, nrow = J-1, ncol = J-1)
    A[upper.tri(A)] <- 0
    A <- kronecker(A, Diagonal(nSam, 1))

    mu.hat <- matrix(A %*% c(pi.hat[, 1: (J-1)]), nrow = nSam)
    y.tilde <- matrix(A %*% c(y.matrix[, 1: (J-1)]), nrow = nSam)

    # W <- get_V_POM.1101(mu.hat, J = J, nSam = nSam)
    W <- get_V_POM(mu.hat, J = J, nSam = nSam)
    D <- quickInverse(W, dimen = J)
    resid <- y.tilde - mu.hat

    resid.list <- list()
    for (i in 1: (J - 1)){
      # resid.list[[paste0('y', i)]] <- c(resid[(1 + (i - 1) * nSam) : (i * nSam)])
      resid.list[[paste0('y', i)]] <- resid[, i]
    }
    y.star.list <- list(unlist(blockdiag.prod(resid.list, D, nSam)))
    W.list <- list(W)
  } else{
    stop('data.type should be either \"nominal\" or \"ordinal\".')
  }

  return(list(X = X, y.star.list = y.star.list, W.list = W.list, J = J, nSam = nSam))
}


#' Sequence Kernel Associaiton Test for Multi-Categorical outcomes (SKAT-MC)
#'
#' @param null.model an output object from the `SKATMC_Null_Model` function.
#' @param G a numeric genotype matrix with each row as a different individual
#' and each column as a separate variant Each genotype should be coded as 0, 1, 2,
#' for AA, Aa, aa, where A is a major allele and a is a minor allele. Missing values
#' are teporarily not allowed.
#' @param weights a numeric vector of weights for the weighted kernels.
#' It is \eqn{\sqrt(w)} in the SKAT paper. We do not recommend to use weights
#' when analyzing common variants.
#' @param weights.beta a numeric vector of parameters for the beta weights for the
#' weighted kernels with default (1, 25). If you want to use your own weights, please use the “weights”
#' parameter. It will be ignored if “weights” parameter is not null.
#'
#' @return
#' Returns a vector of p values with each reference, and if more than one p value
#' generated, the function will automatically return a overall p value from the
#' Cauchy combination. For nominal data. Returns a single p value for ordinal data.
#' @export
#'
#' @examples
SKATMC <- function(null.model, G, weights = NULL, weights.beta = c(1, 10)){
  # browser()
  ###
  # no weights, weights = F
  # beta weights, weights = NULL or T, specify weights.beta
  # other weights, specify weights
  ###

  nSam <- nrow(G)
  nSNP <- ncol(G)
  if (any(is.na(G))) stop('missing values in the genotype matrix.')
  if (nSam != null.model$nSam) stop('Sample size of genotype matrix and
                                          covariates/outcomes do not match.')
  if (!is.null(weights) & !weights) {
    weights <- rep(1, nSNP)
  } else if (is.null(weights) | isTRUE(weights)){
    maf <- colSums(G)/(2 * nSam)
    if (any(maf > 1 | maf < 0)) stop('Minor allele frequency is greater than 1 or
                                     less than 0. Check the genotype matrix.')
    if (any(maf > 0.5)) {
      # warning('Some SNPs have minor allele frequency greater than 0.5, will be flipped to less than 0.5.')
      G[, maf > 0.5] <- 2 - G[, maf > 0.5]
      maf[maf > 0.5] <- 1 - maf[maf > 0.5]
    }
    weights <- dbeta(maf, weights.beta[1], weights.beta[2])
  } else if (length(weights) != nSNP){
    stop('Check the length of weights.')
  }

  if (nSNP > 1000) cat('Number of SNPs is greater than 1000, algorithm might be slow.\n')
  G <- t(G) # by convension, G should be a n * p matrix, but the algorithm needs a p * n matrix
  X <- null.model$X
  y.star.list <- null.model$y.star.list
  W.list <- null.model$W.list
  J <- null.model$J
  all.ref.levels <- length(W.list)
  pv.list <- rep(0, all.ref.levels)
  names(pv.list) <- names(W.list)


  for (i in seq_len(all.ref.levels)){
    tXWX.I <- get_tXWXI(X = X, W = W.list[[i]], J = J)
    qua <- get_qua_linear(w = weights, G = G, W = W.list[[i]], X = X, c = tXWX.I, J = J)
    # lambda <- eigen(qua, symmetric = T, only.values = T) # eigenvalues for chi-square
    score.stat <- get_score.stat(y.star = y.star.list[[i]], W = W.list[[i]], w = weights, G = G, J = J, nSam)
    p.davies <- Get_Davies_PVal(score.stat/2, qua)$p.value
    pv.list[i] <- p.davies
  }

  if (all.ref.levels > 1){
    p.cauchy <- 1 - pcauchy(mean(tan((0.5 - pv.list)*pi)), location = 0, scale = 1)
    pv.list <- c(pv.list, p.omnibus = p.cauchy)
  }
  # browser()

  # p.davies1 <- CompQuadForm::davies(score.stat, lambda$values, lim = 1000, acc = 10^-6)$Qq
  # return(p.davies)
  return(pv.list)
}
