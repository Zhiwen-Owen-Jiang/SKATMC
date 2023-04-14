# These functions are designed for SKAT-MC. Currently, only linear kernel is available.

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









