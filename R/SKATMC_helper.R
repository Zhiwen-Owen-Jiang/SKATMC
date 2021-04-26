# These functions are designed for SKAT-MC. Currently, only linear kernel is available.



getPvalueLOGIT <- function(X, df, y.matrix, w = NULL, G, outfile.logit, J){
  nSam <- dim(df)[1]

  invisible(capture.output(fit <- glm(y.matrix ~ v1 + v2, data = df, family = 'binomial')))
  mu.hat <- fit$fitted.values
  rm(fit)

  V <- Diagonal(nSam, (mu.hat * (1 - mu.hat)))
  y.star <- c(y.matrix[, 1]) - c(mu.hat)
  rm(df)

  ## mixture chi-square distribution and test statistic
  tXV <- crossprod(X, V)
  tXVXI <- solve(crossprod(X, t(tXV)))
  P0 <- V - crossprod(tXV, crossprod(tXVXI, tXV))
  qua <- tcrossprod(tcrossprod(G, P0), G)
  lambda <- eigen(qua, symmetric = T, only.values = T) # eigenvalues for chi-square
  score.stat <- crossprod(crossprod(t(G), y.star))

  p.davies <- davies(score.stat, lambda$values, acc = 10^(-6))$Qq
  rm(lambda)

  if (!is.null(outfile.logit)){
    output(p.davies, outfile.logit)
  } else {
    return(p.davies)
  }

}




getPvalueGLM0920.power <- function(X, y.matrix, df, w, G, outfile.glm, J){
  browser()
  nSam <- dim(df)[1]
  ## fit a generalized logit model
  df$subtype <- factor(df$subtype, levels = as.character(c(J, c(1: (J - 1)))))
  invisible(capture.output(fit <- multinom(subtype ~ ., data = df)))

  ## fitted value and regression coefficients
  mu.hat <- fit$fitted.values[, c(c(2: J), 1)]  ## reference level should be at the last column
  # theta.hat <- c(t(coef(fit)))
  ## D.hat W.hat and y.star
  D <- get_D_GLM.1011(mu.hat, J = J, nSam)
  # X.matrix <- kronecker(diag(1, (J - 1)), X)
  # y.star <- D * (c(y.matrix[, c(1: (J - 1))]) - c(mu.hat[, c(1: (J - 1))])) + X.matrix %*% theta.hat
  y.star <- D * (c(y.matrix[, c(1: (J - 1))]) - c(mu.hat[, c(1: (J - 1))]))

  WI <- get_WI_GLM.1011(mu.hat, D = D, J = J, nSam)
  W <- quickInverse.1011(WI, dimen = J) # should be equal to J

  # W <- get_V_GLM.1101(mu = mu.hat, J = J)
  # W.matrix <- get_matrix_W(W = W, J = J)

  tXWX.I <- get_c(X = X, W = W, J = J)

  ## y.star - Xbeta.hat
  # tXW <- get_tXW(X = X, W = W, J = J)
  # resid <- y.star - X.matrix %*% c %*% (tXW %*% y.star)

  ## mixture chi-square distribution and test statistic
  qua <- get_qua_linear_0211(w = w, G = G, W = W, X = X, c = tXWX.I, J = J)
  lambda <- eigen(qua, symmetric = T, only.values = T) # eigenvalues for chi-square
  score.stat <- get_score.stat_0211(y.star. = y.star, W. = W, w. = w, G. = G, J. = J, nSam)

  p.davies <- CompQuadForm::davies(score.stat, lambda$values, acc = 10^-6)$Qq
  # p.davies <- 1 - sw(coeff = lambda$values, x = score.stat)

  if (!is.null(outfile.glm)){
    output(p.davies, outfile.glm)
  } else {
    return(p.davies)
  }
}



getPvaluePOM1127.power <- function(X, y.matrix, df, w, G, outfile.pom, J){
  nSam <- dim(df)[1]
  ## fit a proportional odds model
  df$subtype <- factor(df$subtype, levels = as.character(c(1: J)))
  fit.pom <- polr(subtype ~ ., data = df, method = 'logistic')
  pi.hat.pom <- fit.pom$fitted.values
  ## regression coefficients
  # theta.hat <- c(fit.pom$zeta, fit.pom$coefficients)

  ## cumulative responses and probabilities
  mu.hat.pom <- c(pi.hat.pom[, 1])
  j = 2
  while (j < J){
    mu.hat.i <- c(pi.hat.pom[, 1])
    for (i in 2:j){
      mu.hat.i <- mu.hat.i + pi.hat.pom[, i]
    }
    mu.hat.pom <- c(mu.hat.pom, mu.hat.i)
    j = j + 1
  }
  mu.hat.pom <- matrix(mu.hat.pom, ncol = J - 1)

  ## cumulative outcomes
  y.tilde.pom <- c(y.matrix[, 1])
  j = 2
  while (j < J){
    y.tilde.i <- c(y.matrix[, 1])
    for (i in 2:j){
      y.tilde.i <- y.tilde.i + y.matrix[, i]
    }
    y.tilde.pom <- c(y.tilde.pom, y.tilde.i)
    j = j + 1
  }

  # X.matrix <- kronecker(diag(1, (J - 1)), X)
  # X.new <- cbind(bdiag(replicate(J-1, rep(1, nSam), simplify = F)),
  #                 rep(-X[, 2], J-1), rep(-X[, 3], J-1))

  W.pom <- get_V_POM.1101(mu.hat.pom, J = J, nSam)
  D.pom <- quickInverse.1011(W.pom, dimen = J, nSam)
  # W.matrix <- get_matrix_W(W = W.pom, J = J)
  # D.pom <- get_matrix_W(W = D.pom, J = J)
  resid.pom <- y.tilde.pom - c(mu.hat.pom)

  resid.pom.list <- list()
  for (i in 1: (J - 1)){
    resid.pom.list[[paste0('y', i)]] <- c(resid.pom[(1 + (i - 1) * nSam) : (i * nSam)])
  }
  y.star.pom <- unlist(blockdiag.prod(resid.pom.list, D.pom, nSam))

  # tXWX.I <- solve(crossprod(X.new, crossprod(W.matrix, X.new)))
  tXWX.I <- get_c(X = X, W = W.pom, J = J)

  ## mixture chi-square and test statistic
  qua.pom <- get_qua_linear_0211(w = w, G = G, W = W.pom, X = X, c = tXWX.I, J = J)
  lambda.pom <- eigen(qua.pom, symmetric = T, only.values = T)
  score.stat.pom <- get_score.stat_0211(y.star = y.star.pom, W = W.pom, w = w, G = G, J = J, nSam)

  ## p value
  p.davies.pom <- davies(score.stat.pom, lambda.pom$values, acc = 10^-6)$Qq

  if (!is.null(outfile.pom)){
    output(p.davies.pom, outfile.pom)
  } else {
    return(p.davies.pom)
  }
}



getPvalueburden <- function(X, df, b, G, outfile.burden, J, model){
  if (model == 'glm'){
    # df$y.model <- factor(df$y.model, levels = as.character(c(J, c(1: (J - 1)))))
    invisible(capture.output(fit <- multinom(y.model ~ v1 + v2 + t(G) %*% b, data = df)))
    # print(fit)
    # return(fit)
    if (J == 3){
      test.fit <- linearHypothesis(fit, hypothesis.matrix = matrix(c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1), nrow = 2, byrow = T))
    } else if (J == 5){
      test.fit <- linearHypothesis(fit, hypothesis.matrix = matrix(c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1), nrow = 4, byrow = T))
    }

    if (!is.null(outfile.burden)){
      output(test.fit$`Pr(>Chisq)`[2], outfile.burden)
    } else {
      return(test.fit$`Pr(>Chisq)`[2])
    }
  }
  if (model == 'pom'){
    df$y.model <- factor(df$y.model, levels = as.character(c(1: J)))
    fit.pom <- polr(y.model ~ v1 + v2 + t(G) %*% b, data = df, method = 'logistic', Hess = T)
    test.fit.pom <- linearHypothesis(fit.pom, hypothesis.matrix = matrix(c(0,0,1), nrow = 1, byrow = T))
    if (!is.null(outfile.burden)){
      output(test.fit.pom$`Pr(>Chisq)`[2], outfile.burden)
    } else {
      return(test.fit.pom$`Pr(>Chisq)`[2])
    }
  }
}


getPvalueburden.score <- function(X, df, G, outfile.burden, J, nSam){
  # df$y.model <- factor(df$y.model, levels = as.character(c(J, c(1: (J - 1)))))
  invisible(capture.output(fit.burden <- multinom(y.model ~ v1 + v2, data = df)))
  X.k <- kronecker(diag(J-1), X)
  X2 <- colSums(G)
  X2.k <- kronecker(diag(J-1), X2)
  y.matrix <- model.matrix(~ 0 + y.model, data = df)
  mu.hat <- fit.burden$fitted.values
  # score.fun <- t(X) %*% (c(y.matrix[, 2: 5] - mu.hat[, 2: 5]))
  V <- get_WI_GLM.1011(mu = mu.hat[, c(2:J, 1)], D = rep(1, nSam * (J-1)), J = J, nSam = nSam)
  tXVX.I <- get_c(X, W = V, J = J)
  tXVX2 <- get_c2(X, W = V, X2, J = J)
  tX2VX2 <- solve(get_c(X2, W = V, J = J))
  E2 <- X2.k - X.k %*% tXVX.I %*% tXVX2
  # score.fun <- t(E2) %*% (c(y.matrix[, 2: J] - mu.hat[, 2: J]))
  score.fun <- t(E2) %*% (c(y.matrix[, 2:J] - mu.hat[, 2: J]))
  info <- tX2VX2 - t(tXVX2) %*% tXVX.I %*% tXVX2
  chisq.stat <- t(score.fun) %*% solve(info) %*% score.fun
  print(chisq.stat)
  pv <- pchisq(chisq.stat, J-1, lower.tail = F)
  if (!is.null(outfile.burden)){
    output(pv, outfile.burden)
  } else {
    return(pv)
  }
}




blockdiag.prod <- function(A, B, nSam = nSam){
  ## A and B should be block diagonal matrix, stored in form of list containing
  ## only upper triangle (because here we are dealing with symmetric matrix)
  ## A should be a row vector, B should be square matrix
  out <- list()
  lenA <- length(A)
  # lenB <- length(B)
  for (i in 1:lenA){
    out[[paste0('w', i)]] <- rep(0, nSam)
    for (j in 1:lenA){
      if (i == 1 | i == 2){
        out[[paste0('w', i)]] <- out[[paste0('w', i)]] + A[[j]] * B[[i + sum(0:(j-1))]]
      } else{
        if (j <= i){
          out[[paste0('w', i)]] <- out[[paste0('w', i)]] + A[[j]] * B[[sum(0:(i-1)) + j]]
        } else{
          out[[paste0('w', i)]] <- out[[paste0('w', i)]] + A[[j]] * B[[sum(0:(i-1)) + i + sum(i:(j - 1))]]
        }
      }
    }
  }
  return(out)
}


matrixlist.prod <- function(A, B, type, nSam = nSam){
  ## product of two 'vectors' in form of list
  ## type 1 represents a row (1*n) times a column (n*1)
  ## type 2 represents a column (n*1) times a row (1*n)
  ## type 3 represents a number (1*1) times a row (1*n) or a column (n*1) times
  ## a number (1*1)
  lenA <- length(A)
  lenB <- length(B)
  if (lenA == lenB & type == 1){
    out <- rep(0, nSam)
    for (i in 1:lenA){
      out <- out + A[[i]] * B[[i]]
    }
    return(list(w1 = out))
  } else if (lenA == lenB & type == 2){
    out <- list()
    for (i in 1:lenA){
      for (j in 1:i){
        if (i == j){
          out[[paste0('w', i, j)]] <- A[[i]] * B[[j]]
        } else{
          out[[paste0('w', j, i)]] <- A[[i]] * B[[j]]
        }
      }
    }
    return(out)
  }


  if (lenA == 1 & type == 3){
    out <- list()
    for (i in 1:lenB){
      out[[paste0('w', i, lenB)]] <- A[[1]] * B[[i]]
    }
    return(out)
  }

  if (lenB == 1 & type == 3){
    out <- list()
    for (i in 1:lenA){
      out[[paste0('w', i, lenA)]] <- A[[i]] * B[[1]]
    }
    return(out)
  }
}



quickInverse <- function(M, dimen, nSam = nSam){
  ## suppose M is a square matrix
  length.M <- length(M)

  if (length.M != 3){
    Ai <- quickInverse(M[c(1: sum(1:(dimen - 2)))], dimen - 1)
  }else {
    Ai <- list(w11 = 1/unlist(M[1]))
  }

  B <- M[c(sum(1: (dimen - 2)) + 1): c(sum(1: (dimen - 2)) + dimen - 2)]
  D <- unlist(M[length.M])

  BAi <- blockdiag.prod(B, Ai, nSam)
  DI <- list()
  DI[[paste0('w', dimen, dimen)]] <- 1/(D - unlist(matrixlist.prod(BAi, B, type = 1)))
  BI <- lapply(matrixlist.prod(BAi, DI, type = 3), function(x) -1 * x)
  AI <- mapply('-', Ai, matrixlist.prod(BI, BAi, type = 2), SIMPLIFY = F)
  MI <- c(AI, BI, DI)
  return(MI)
}


quickInverse.1011 <- function(M, dimen, nSam = nSam){
  ## suppose M is a square matrix
  length.M <- length(M)

  if (length.M != 3){
    Ai <- quickInverse.1011(M[c(1: sum(1:(dimen - 2)))], dimen - 1)
  }else {
    Ai <- list(w11 = 1/unlist(M[1]))
  }

  B <- M[c(sum(1: (dimen - 2)) + 1): c(sum(1: (dimen - 2)) + dimen - 2)]
  D <- unlist(M[length.M])

  BAi <- blockdiag.prod(B, Ai, nSam)
  DI <- list()
  DI[[paste0('w', dimen, dimen)]] <- 1/(D - unlist(matrixlist.prod(BAi, B, type = 1)))
  BI <- lapply(matrixlist.prod(BAi, DI, type = 3), function(x) -1 * x)
  AI <- mapply('-', Ai, matrixlist.prod(BI, BAi, type = 2), SIMPLIFY = F)
  MI <- c(AI, BI, DI)
  return(MI)
}

get_WI_GLM.1011 <- function(mu, D, J, nSam = nSam){
  ## reference level should be put in the last column
  vv <- rep(0, (J - 1) * nSam)
  vc <- rep(0, sum(1: (J - 2)) * nSam)
  n <- 1
  m <- 1
  for (j in c(1:(J-1))){
    for (k in c(1:(J-1))){
      for (i in c(1:nSam)) {
        if (j == k){
          num <- mu[i,j]*(1-mu[i,j])
          vv[n] <- num
          n = n + 1
        }
        else if (j > k){
          num <- -mu[i,j]*mu[i,k]
          # num <- 0
          vc[m] <- num
          m = m + 1
        }
      }
    }
  }

  wiv <- D^2 * vv

  j <- 1
  cc <- c()
  dd <- c()
  # while (j < J-1) {
  #   dd <- c(dd, rep(D[((j - 1) * nSam + 1): (j * nSam)], J-1-j))
  #   cc <- c(cc, D[-(1: (j * nSam))])
  #   j = j + 1
  # }

  while (j < J-1) {
    dd <- c(dd, D[1: (j * nSam)])
    cc <- c(cc, rep(D[(j * nSam + 1): ((j + 1) * nSam)], j))
    j = j + 1
  }

  wic <- dd * vc * cc

  WI <- list()
  k = 1
  for (i in c(1: (J-1))){
    for (j in c(1 : i)){
      if (i == j) {
        WI[[paste0('w',i,j)]] = wiv[((j - 1) * nSam + 1): (j * nSam)]
      }else{
        WI[[paste0('w',j,i)]] = wic[((k - 1) * nSam + 1): (k * nSam)]
        k = k + 1
      }
    }
  }

  return(WI)
}

get_D_GLM.1011 <- function(mu, J, nSam = nSam){
  ## reference level should be put in the last column
  d <- rep(0, (J-1)*nSam)
  n <- 1
  for (j in c(1:(J-1))){
    for (i in c(1:nSam)) {
      # num <- 1/(mu[i,j]*(1-mu[i,j]))
      # num <- 1/(mu[i,j])
      num <- (mu[i,j] + mu[i,J])/(mu[i,j]*mu[i,J])
      # num <- 1
      d[n] <- num
      n = n + 1
    }
  }
  return(d)
}


get_WI_GLM <- function(mu, D, J, nSam = nSam){
  # browser()
  vc <- matrix(0, nrow = nSam, ncol = sum(1: (J - 2)))
  m = 1
  vv <- mu * (1 - mu)

  for (j in c(1:(J-1))){
    for (k in c(1:(J-1))){
        if (j > k){
          vc[, m] <- -mu[, j] * mu[, k]
          m = m + 1
      }
    }
  }

  wiv <- matrix(D^2 * c(vv), nrow = nSam)
  D.m <- matrix(D, nrow = nSam)

  dd <- c(D.m[, unlist(mapply(seq, 1, 1: (J - 2)))])
  cc <- c(D.m[, unlist(mapply(rep, 2: (J - 1), 1: (J - 2)))])

  wic <- matrix(dd * vc * cc, nrow = nSam)

  WI <- list()
  k = 1
  for (i in c(1: (J-1))){
    for (j in c(1 : i)){
      if (i == j) {
        WI[[paste0('w',i,j)]] = wiv[, j]
      }else{
        WI[[paste0('w',j,i)]] = wic[, k]
        k = k + 1
      }
    }
  }

  return(WI)
}



get_V_POM.1101 <- function(mu, J, nSam = nSam){
  ## reference level should be put in the last column
  vv <- rep(0, (J - 1) * nSam)
  #  <- rep(1, (J - 1) * nSam)
  vc <- rep(0, sum(1: (J - 2)) * nSam)
  n <- 1
  m <- 1

  for (j in c(1:(J-1))){
    for (k in c(1:(J-1))){
      for (i in c(1:nSam)) {
        if (j == k){
          num <- mu[i,j]*(1-mu[i,j])
          vv[n] <- num
          n = n + 1
        }
        else if (j > k){
          num <- mu[i,k]*(1 - mu[i,j])
          vc[m] <- num
          m = m + 1
        }
      }
    }
  }

  # wiv <- vv
  # j <- 1
  # cc <- c()
  # dd <- c()
  # while (j < J-1) {
  #   dd <- c(dd, rep(D[((j - 1) * nSam + 1): (j * nSam)], J-1-j))
  #   cc <- c(cc, D[-(1: (j * nSam))])
  #   j = j + 1
  # }

  # while (j < J-1) {
  #   dd <- c(dd, I[1: (j * nSam)])
  #   cc <- c(cc, rep(I[(j * nSam + 1): ((j + 1) * nSam)], j))
  #   j = j + 1
  # }

  # wic <- dd * vc * cc
  # wic <- vc
  V <- list()
  k = 1
  for (i in c(1: (J-1))){
    for (j in c(1 : i)){
      if (i == j) {
        V[[paste0('w',i,j)]] = vv[((j - 1) * nSam + 1): (j * nSam)]
      }else{
        V[[paste0('w',j,i)]] = vc[((k - 1) * nSam + 1): (k * nSam)]
        k = k + 1
      }
    }
  }

  return(V)
}


get_V_POM <- function(mu, J, nSam){
  # browser()
  vc <- matrix(0, nrow = nSam, ncol = sum(1: (J - 2)))
  m = 1
  vv <- mu * (1 - mu)

  for (j in c(1:(J-1))){
    for (k in c(1:(J-1))){
      if (j > k){
        vc[, m] <- mu[, k] * (1 - mu[, j])
        m = m + 1
      }
    }
  }

  V <- list()
  k = 1
  for (i in c(1: (J-1))){
    for (j in c(1 : i)){
      if (i == j) {
        V[[paste0('w',i,j)]] = vv[, j]
      } else{
        V[[paste0('w',j,i)]] = vc[, k]
        k = k + 1
      }
    }
  }


  return(V)
}


get_matrix_W <- function(W, J){
  W.matrix <- c()
  for (i in 1: (J - 1)){
    W.i <- diag(W[[1 + sum(0:(i-1))]])
    for (j in 2: (J - 1)){
      if (i == 1 | i == 2){
        W.i <- cbind(W.i, diag(W[[i + sum(0:(j-1))]]))
      } else {
        if (j <= i){
          W.i <- cbind(W.i, diag(W[[sum(0:(i-1)) + j]]))
        } else {
          W.i <- cbind(W.i, diag(W[[sum(0:(i-1)) + i + sum(i:(j - 1))]]))
        }
      }
    }
    W.matrix <- rbind(W.matrix, W.i)
  }
  return(W.matrix)
}


get_tXWXI <- function(X, W, J){
  tXWX.list <- lapply(W, function(w){crossprod(w * X, X)})
  c.matrix <- c()
  for (i in 1: (J - 1)){
    c.matrix.i <- tXWX.list[[1 + sum(0:(i-1))]]
    for (j in 2: (J - 1)){
      if (i == 1 | i == 2){
        c.matrix.i <- cbind(c.matrix.i, tXWX.list[[i + sum(0:(j-1))]])
      } else {
        if (j <= i){
          c.matrix.i <- cbind(c.matrix.i, tXWX.list[[sum(0:(i-1)) + j]])
        } else {
          c.matrix.i <- cbind(c.matrix.i, tXWX.list[[sum(0:(i-1)) + i + sum(i:(j - 1))]])
        }
      }
    }
    c.matrix <- rbind(c.matrix, c.matrix.i)
  }
  return(solve(c.matrix))
}



get_c <- function(X, W, J){
  tXWX.list <- lapply(W, function(w){crossprod(w * X, X)})
  c.matrix <- c()
  for (i in 1: (J - 1)){
    c.matrix.i <- tXWX.list[[1 + sum(0:(i-1))]]
    for (j in 2: (J - 1)){
      if (i == 1 | i == 2){
        c.matrix.i <- cbind(c.matrix.i, tXWX.list[[i + sum(0:(j-1))]])
      } else {
        if (j <= i){
          c.matrix.i <- cbind(c.matrix.i, tXWX.list[[sum(0:(i-1)) + j]])
        } else {
          c.matrix.i <- cbind(c.matrix.i, tXWX.list[[sum(0:(i-1)) + i + sum(i:(j - 1))]])
        }
      }
    }
    c.matrix <- rbind(c.matrix, c.matrix.i)
  }
  return(solve(c.matrix))
}


get_c2 <- function(X1, W, X2, J){
  tXWX.list <- lapply(W, function(w){crossprod(w * X1, X2)})
  c.matrix <- c()
  for (i in 1: (J - 1)){
    c.matrix.i <- tXWX.list[[1 + sum(0:(i-1))]]
    for (j in 2: (J - 1)){
      if (i == 1 | i == 2){
        c.matrix.i <- cbind(c.matrix.i, tXWX.list[[i + sum(0:(j-1))]])
      } else {
        if (j <= i){
          c.matrix.i <- cbind(c.matrix.i, tXWX.list[[sum(0:(i-1)) + j]])
        } else {
          c.matrix.i <- cbind(c.matrix.i, tXWX.list[[sum(0:(i-1)) + i + sum(i:(j - 1))]])
        }
      }
    }
    c.matrix <- rbind(c.matrix, c.matrix.i)
  }
  return(c.matrix)
}


get_matrix_wsqrtGW <- function(wsqrtGW.list, J){
  wsqrtGW <- c()
  for (i in 1: (J - 1)){
    wsqrtGW.i <- wsqrtGW.list[[1 + sum(0:(i-1))]]
    for (j in 2: (J - 1)){
      if (i == 1 | i == 2){
        wsqrtGW.i <- cbind(wsqrtGW.i, wsqrtGW.list[[i + sum(0:(j-1))]])
      } else {
        if (j <= i){
          wsqrtGW.i <- cbind(wsqrtGW.i, wsqrtGW.list[[sum(0:(i-1)) + j]])
        } else {
          wsqrtGW.i <- cbind(wsqrtGW.i, wsqrtGW.list[[sum(0:(i-1)) + i + sum(i:(j - 1))]])
        }
      }
    }
    wsqrtGW <- rbind(wsqrtGW, wsqrtGW.i)
  }
  return(wsqrtGW)
}



# get_qua_linear <- function(w, G, W, X, c, J){
#   wsqrttG <- sqrt(w) * G
#   wsqrtGW.list <- lapply(W, function(w){t(w * t(wsqrttG))})
#
#   wsqrtGW <- c()
#   for (i in 1: (J - 1)){
#     wsqrtGW.i <- wsqrtGW.list[[1 + sum(0:(i-1))]]
#     for (j in 2: (J - 1)){
#       if (i == 1 | i == 2){
#         wsqrtGW.i <- cbind(wsqrtGW.i, wsqrtGW.list[[i + sum(0:(j-1))]])
#       } else {
#         if (j <= i){
#           wsqrtGW.i <- cbind(wsqrtGW.i, wsqrtGW.list[[sum(0:(i-1)) + j]])
#         } else {
#           wsqrtGW.i <- cbind(wsqrtGW.i, wsqrtGW.list[[sum(0:(i-1)) + i + sum(i:(j - 1))]])
#         }
#       }
#     }
#     wsqrtGW <- rbind(wsqrtGW, wsqrtGW.i)
#   }
#
#   # X.matrix <- kronecker(diag(1, (J - 1)), X)
#   # qua1 <- wsqrtGW %*% kronecker(diag(1, (J - 1)), t(wsqrttG)) - wsqrtGW %*% X.matrix %*% c %*% t(X.matrix) %*% t(wsqrtGW)
#   wsqrtGWX <- wsqrtGW %*% X
#   qua1 <- wsqrtGW %*% kronecker(diag(1, (J - 1)), t(wsqrttG)) - wsqrtGWX %*% c %*% t(wsqrtGWX)
#   return(qua1)
# }



get_qua_linear <- function(w, G, W, X, c, J){
  # browser()
  # weight is sqrt(w)
  twsqrttG <- t(w * G)
  wsqrtGW.list <- lapply(W, function(w1){t(w1 * twsqrttG)})
  wsqrtGWGwsqrt.list <- lapply(wsqrtGW.list, function(x){x %*% twsqrttG})
  wsqrtGWX.list <- lapply(wsqrtGW.list, function(x){x %*% X})
  wsqrtGWGwsqrt <- get_matrix_wsqrtGW(wsqrtGWGwsqrt.list, J)
  wsqrtGWX <- get_matrix_wsqrtGW(wsqrtGWX.list, J)

  qua1 <- wsqrtGWGwsqrt - wsqrtGWX %*% c %*% t(wsqrtGWX)
  return(qua1)
}



get_qua_linear_0211 <- function(w, G, W, X, c, J){
  # browser()
  twsqrttG <- t(sqrt(w) * G)
  wsqrtGW.list <- lapply(W, function(w1){t(w1 * twsqrttG)})
  wsqrtGWGwsqrt.list <- lapply(wsqrtGW.list, function(x){x %*% twsqrttG})
  wsqrtGWX.list <- lapply(wsqrtGW.list, function(x){x %*% X})
  wsqrtGWGwsqrt <- get_matrix_wsqrtGW(wsqrtGWGwsqrt.list, J)
  wsqrtGWX <- get_matrix_wsqrtGW(wsqrtGWX.list, J)

  qua1 <- wsqrtGWGwsqrt - wsqrtGWX %*% c %*% t(wsqrtGWX)
  return(qua1)
}


get_score.stat <- function(y.star, W, w, G, J, nSam = nSam){
  # browser()
  y.star <- matrix(y.star, nrow = nSam)
  y.star.list <- list()
  for (i in 1: (J - 1)){
    # y.star.list[[paste0('y', i)]] <- c(y.star.[(1 + (i - 1) * nSam.) : (i * nSam.)])
    y.star.list[[paste0('y', i)]] <- y.star[, i]
  }
  y.star.W <- blockdiag.prod(y.star.list, W)

  w <- w^2
  twGG <- t(w * G) %*% G
  y.star.WK <- lapply(y.star.W, function(x) {t(x) %*% twGG})
  score.stat <- 0
  for (i in 1:length(y.star.WK)){
    score.stat <- score.stat + y.star.WK[[i]] %*% y.star.W[[i]]
  }
  return(c(score.stat))
}


get_score.stat_0211 <- function(y.star., W., w., G., J., nSam. = nSam){
  # browser()
  y.star.list <- list()
  for (i in 1: (J. - 1)){
    y.star.list[[paste0('y', i)]] <- c(y.star.[(1 + (i - 1) * nSam.) : (i * nSam.)])
  }
  y.star.W <- blockdiag.prod(y.star.list, W.)
  twsqrttG <- t(sqrt(w.) * G.)

  y.star.WG <- lapply(y.star.W, function(x) {t(x) %*% twsqrttG})
  score.stat <- 0
  for (i in 1:length(y.star.WG)){
    score.stat <- score.stat + sum(y.star.WG[[i]]^2)
  }
  return(c(score.stat))
}


output <- function(value, filename){
  sink(filename, append = T)
  cat(paste0(value,'\n'))
  sink()
}



getPvalueGLM0920.tp <- function(X, y.matrix, df, J){
  nSam = dim(df)[1]
  ## fit a generalized logit model
  df$y.model <- factor(df$y.model, levels = as.character(c(J, c(1: (J - 1)))))
  invisible(capture.output(fit <- multinom(y.model ~ ., data = df)))

  ## fitted value and regression coefficients
  mu.hat <- fit$fitted.values[, c(c(2: J), 1)]  ## reference level should be at the last column
  # theta.hat <- c(t(coef(fit)))

  ## D.hat W.hat and y.star
  D <- get_D_GLM.1011(mu.hat, J = J, nSam = nSam)
  X.matrix <- kronecker(diag(1, (J - 1)), X)
  # y.star <- D * (c(y.matrix[, c(1: (J - 1))]) - c(mu.hat[, c(1: (J - 1))])) + X.matrix %*% theta.hat
  y.star <- D * (c(y.matrix[, c(1: (J - 1))]) - c(mu.hat[, c(1: (J - 1))]))

  WI <- get_WI_GLM.1011(mu.hat, D = D, J = J, nSam = nSam)
  W <- quickInverse.1011(WI, dimen = J) # should be equal to J

  # W <- get_V_GLM.1101(mu = mu.hat, J = J)
  # W.matrix <- get_matrix_W(W = W, J = J)
  return(list(X = X, W = W, y.star = y.star))
}



getPvaluePOM1127.tp <- function(X, y.matrix, df, J){
  nSam = dim(df)[1]
  ## fit a proportional odds model
  df$y.model <- factor(df$y.model, levels = as.character(c(1: J)))
  fit.pom <- polr(y.model ~ ., data = df, method = 'logistic')
  pi.hat.pom <- fit.pom$fitted.values
  ## regression coefficients
  # theta.hat <- c(fit.pom$zeta, fit.pom$coefficients)

  ## cumulative responses and probabilities
  mu.hat.pom <- c(pi.hat.pom[, 1])
  j = 2
  while (j < J){
    mu.hat.i <- c(pi.hat.pom[, 1])
    for (i in 2:j){
      mu.hat.i <- mu.hat.i + pi.hat.pom[, i]
    }
    mu.hat.pom <- c(mu.hat.pom, mu.hat.i)
    j = j + 1
  }
  mu.hat.pom <- matrix(mu.hat.pom, ncol = J - 1)

  ## cumulative outcomes
  y.tilde.pom <- c(y.matrix[, 1])
  j = 2
  while (j < J){
    y.tilde.i <- c(y.matrix[, 1])
    for (i in 2:j){
      y.tilde.i <- y.tilde.i + y.matrix[, i]
    }
    y.tilde.pom <- c(y.tilde.pom, y.tilde.i)
    j = j + 1
  }

  X.matrix <- kronecker(diag(1, (J - 1)), X)
  # X.new <- cbind(bdiag(replicate(J-1, rep(1, nSam), simplify = F)),
  #                 rep(-X[, 2], J-1), rep(-X[, 3], J-1))

  W.pom <- get_V_POM.1101(mu.hat.pom, J = J, nSam = nSam)
  D.pom <- quickInverse.1011(W.pom, dimen = J)
  # W.matrix <- get_matrix_W(W = W.pom, J = J)
  # D.pom <- get_matrix_W(W = D.pom, J = J)
  resid.pom <- y.tilde.pom - c(mu.hat.pom)

  resid.pom.list <- list()
  for (i in 1: (J - 1)){
    resid.pom.list[[paste0('y', i)]] <- c(resid.pom[(1 + (i - 1) * nSam) : (i * nSam)])
  }
  y.star.pom <- unlist(blockdiag.prod(resid.pom.list, D.pom))

  return(list(X = X, W = W.pom, y.star = y.star.pom))
}


get.pvalue.gene.GLMPOM <- function(X, W, y.star, w, G, outfile, J){
  # browser()
  nSam = dim(G)[2]
  tXWX.I <- get_c(X = X, W = W, J = J)
  # X.matrix <- kronecker(diag(1, (J - 1)), X)
  ## y.star - Xbeta.hat
  # tXW <- get_tXW(X = X, W = W, J = J)
  # resid <- y.star - X.matrix %*% c %*% (tXW %*% y.star)

  ## mixture chi-square distribution and test statistic
  qua <- get_qua_linear_0211(w = w, G = G, W = W, X = X, c = tXWX.I, J = J)
  lambda <- eigen(qua, symmetric = T, only.values = T) # eigenvalues for chi-square
  score.stat <- get_score.stat_0211(y.star = y.star, W = W, w = w, G = G, J = J, nSam)

  p.davies <- davies(score.stat, lambda$values, lim = 1000, acc = 10^-6)$Qq
  # p.davies <- 1 - sw(lambda$values, score.stat)
  if (!is.null(outfile)){
    output(p.davies, outfile)
  } else {
    return(p.davies)
  }
}




get.pvalue.gene.GLMPOM.DKAT <- function(X, W, y.star, w, G, outfile, J){
  browser()
  nSam = dim(G)[2]
  y.star.list <- list()
  for (i in 1: (J - 1)){
    y.star.list[[paste0('y', i)]] <- c(y.star[(1 + (i - 1) * nSam) : (i * nSam)])
  }
  y.star.W <- blockdiag.prod(y.star.list, W)
  y.star.W.matrix <- do.call(cbind, y.star.W)
  KY <- tcrossprod(y.star.W.matrix)

  sqrtwG <- sqrt(w) * G
  K <- crossprod(sqrtwG)

  p.DKAT <- DKAT(KY, K)
  # p.davies <- 1 - sw(lambda$values, score.stat)
  if (!is.null(outfile)){
    output(p.DKAT, outfile)
  } else {
    return(p.DKAT)
  }
}





getPvalueGLM0920.realdata <- function(X, y.matrix, df, J){
  nSam = dim(df)[1]
  ## fit a generalized logit model
  df$ER <- factor(df$ER, levels = as.character(c(J-1, c(0: (J - 2)))))
  invisible(capture.output(fit <- multinom(ER ~ ., data = df)))

  ## fitted value and regression coefficients
  mu.hat <- fit$fitted.values[, c(c(2: J), 1)]  ## reference level should be at the last column
  # theta.hat <- c(t(coef(fit)))

  ## D.hat W.hat and y.star
  D <- get_D_GLM.1011(mu.hat, J = J, nSam = nSam)
  # X.matrix <- kronecker(diag(1, (J - 1)), X)
  # y.star <- D * (c(y.matrix[, c(1: (J - 1))]) - c(mu.hat[, c(1: (J - 1))])) + X.matrix %*% theta.hat
  y.star <- D * (c(y.matrix[, c(1: (J - 1))]) - c(mu.hat[, c(1: (J - 1))]))

  WI <- get_WI_GLM.1011(mu.hat, D = D, J = J, nSam = nSam)
  # W <- quickInverse.1011(WI, dimen = J) # should be equal to J

  # W <- get_V_GLM.1101(mu = mu.hat, J = J)
  # W.matrix <- get_matrix_W(W = W, J = J)
  return(list(X = X, WI = WI, y.star = y.star))
}







