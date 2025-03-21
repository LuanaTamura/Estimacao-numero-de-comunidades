######
#Teste de hipóteses
######

#vetor de probabilidade q=pi
#n número de vértices 
sample.z <- function(n,q){
  k <- length(q)
  return(sample(x = seq(1:k), n, replace = T, prob = q))
}

#P será a matriz B
sample.g <- function(z,P){
  n <- length(z)
  mat.adj <- matrix(0,n,n)
  k <- dim(P)[1]
  for(i in 1:n){
    for(j in i:n){
      mat.adj[i,j] <- mat.adj[j,i] <- rbinom(1,1,P[z[i],z[j]])
    }
  }
  diag(mat.adj) <- 0
  return(mat.adj)
}

test.hip <- function(n){
  
  Generate.A <- function(P) {
    # generates random adjacency matrix P
    n <- dim(P)[1]
    upper.tri.ind <- upper.tri(P)
    p.upper <- P[upper.tri.ind]
    A.upper <- rbinom(n*(n-1)/2, 1, p.upper)
    A <- matrix(0, ncol = n, nrow = n)
    A[upper.tri.ind] <- A.upper
    A <- A + t(A)
    return(A)
  }
  
  #estimar as comunidades z
  SpecClust <- function(A, K, sphere = F, n.dim = K, lift = 0, nstart = 20) {
    # spectral clustering
    if (is.null(n.dim)) {n.dim <- K}
    n <- nrow(A)
    A <- A + matrix(lift, ncol = n, nrow = n)
    if (K == 1) {
      return(rep(1, n))
    }
    U <- slanczos(A, n.dim)$vectors
    if (sphere) {
      for (i in 1:n) {
        U[i,] = U[i,] / sqrt(sum(U[i,]^2))
      }
    }
    return(kmeans(U, K, nstart = nstart)$cluster)
    #  return(clusters(kcca(U, K)))
  }
  
  #vetor z que dará as comunidades de cada vértice
  ClustVec2Mat <- function(clust.vec, K) {
    # convert a membership vector to a membership matrix, the number of different values in
    # clust.vec must be at most K
    clust.mat <- matrix(0, ncol = K, nrow = length(clust.vec))
    for (i in 1:K) {
      clust.mat[clust.vec==i,i] <- 1
    }
    return(clust.mat)
  }
  
  #quanto de erro cometido do estimado e o verdadeiro (z^ e o vdd)
  ClustError <- function(T1, T2) {
    # cluster error calculation, taking into account of possibly label permutation
    # Inputs:
    #   T1 - a *vector* of estimated cluster membership
    #   T2 - the true membership *matrix*
    n <- dim(T2)[1]
    K <- dim(T2)[2]
    T1.new <- ClustVec2Mat(T1, K)
    perm.list <- permn(1:K)
    l <- length(perm.list)
    rec.err <- rep(0, l)
    for (i in 1:l) {
      rec.err[i] <- sum(sum(abs(T1.new[, perm.list[[i]] ] - T2))) / 2
    }
    i.star <- which(rec.err==min(rec.err))[1]
    return(list(clust.est = T1.new %*% diag(K)[perm.list[[i.star]],], err = rec.err[i.star]))
  }
  
  #estimar o B
  EmpB <- function(A, clusters) {
    # empirical connectivity matrix
    values <- sort(unique(clusters))
    K <- length(values)
    if (K == 1) {
      n <- nrow(A)
      B = sum(A) / (n * (n - 1))
    } else {
      B <- matrix(0, nrow = K, ncol = K)
      clust.ind <- list()
      clust.size <- rep(0, K)
      for (i in 1:K) {
        clust.ind[[i]] <- which(clusters == i)
        clust.size[i] <- length(clust.ind[[i]])
      }
      for (i in 1:(K-1)) {
        B[i,i] <- ifelse(clust.size[i] < 2, 0.5,
                         (sum(A[clust.ind[[i]], clust.ind[[i]]]) + 2) / 
                           (clust.size[i]^2 - clust.size[i] + 2))
        for (j in (i+1):K) {
          B[i,j] <- (sum(A[clust.ind[[i]], clust.ind[[j]]]) + 1) / 
            (clust.size[i] * clust.size[j] + 1)
          B[j,i] <- B[i,j]
        }
      }
      B[K,K] <- ifelse(clust.size[K] < 2, 0.5,
                       (sum(A[clust.ind[[K]], clust.ind[[K]]]) + 2) /
                         (clust.size[K]^2 - clust.size[K] + 2))
    }
    return(B)
  }
  
  #estatística do teste 
  GoFStat <- function(A, K0) {
    n <- nrow(A)
    clusters <- SpecClust(A, K0) #z^
    B.hat <- EmpB(A, clusters) 
    Theta.hat <- Generate.theta(n, K0, table(clusters)) #matriz para auxílio nos calculos
    P.hat <- Theta.hat %*% B.hat %*% t(Theta.hat)
    tilde.A <- (A - P.hat + diag(diag(P.hat))) / sqrt((n-1)*P.hat*(1-P.hat))
    lambda <- slanczos(tilde.A, k = 1, kl = 1)$values
    lambda.1 <- lambda[1]
    lambda.n <- -lambda[2]
    sig.1 <- max(abs(lambda))
    return(n^(2/3) * (c(lambda.1, lambda.n, sig.1) - 2))
  }
  
  GoFTestBoot <- function(A, K0 = 1, n.dim = K0, n.b = 50, use.boot = T) {
    n <- nrow(A)
    clusters <- SpecClust(A, K0, n.dim = n.dim)
    B.hat <- EmpB(A, clusters)
    Theta.hat <- ClustVec2Mat(clusters, K0)
    P.hat <- Theta.hat %*% B.hat %*% t(Theta.hat)
    tilde.A <- (A - P.hat + diag(diag(P.hat))) / sqrt((n-1)*P.hat*(1-P.hat))
    lambda <- slanczos(tilde.A, k = 1, kl = 1)$values
    lambda.1 <- lambda[1]
    lambda.n <- -lambda[2]
    test.stat.boot <- NULL
    p.value.boot <- NULL
    if (use.boot) {
      mu.tw <- -1.206533574582
      sd.tw <- sqrt(1.607781034581)
      lambda.1.vec <- rep(0, n.b)
      lambda.n.vec <- rep(0, n.b)
      for (i in 1:50) {
        A.i <- Generate.A(P.hat)
        tilde.A.i <- (A.i - P.hat + diag(diag(P.hat))) /
          sqrt((n - 1) * P.hat * (1 - P.hat))
        lambda.i <- slanczos(tilde.A.i, k = 1, kl = 1)$values
        lambda.1.vec[i] <- lambda.i[1]
        lambda.n.vec[i] <- -lambda.i[2]
      }
      mu.lambda.1.hat <- mean(lambda.1.vec, na.rm = T)
      sig.lambda.1.hat <- sd(lambda.1.vec, na.rm = T)
      lambda.1.prime <- mu.tw + (lambda.1 - mu.lambda.1.hat) /
        sig.lambda.1.hat * sd.tw
      mu.lambda.n.hat <- mean(lambda.n.vec, na.rm = T)
      sig.lambda.n.hat <- sd(lambda.n.vec, na.rm = T)
      lambda.n.prime <- mu.tw + (lambda.n - mu.lambda.n.hat) /
        sig.lambda.n.hat * sd.tw
      test.stat.boot <- c(lambda.1.prime, lambda.n.prime)
    }
    return(list(test.stat = c(n^(2/3) * (lambda.1 - 2),
                              n^(2/3) * (lambda.n - 2)),
                test.stat.boot = test.stat.boot))
  }
  
  SeqTestSBM <- function(A, alpha = 0.0001, K.min = 1, K.max = 5, n.dim = NULL, n.b = 50, use.boot = T) {
    n <- nrow(A)
    threshold <- qtw(alpha / 2, lower.tail = F) #quantil da tracy widom
    n.alpha <- length(alpha)
    K.est <- rep(NA, n.alpha)
    K.est.boot <- K.est
    if (use.boot) {
      for (K0 in K.min:K.max) {
        test.stat.K0 <- GoFTestBoot(A, K0 = K0, n.dim = n.dim, n.b = n.b, use.boot = use.boot)
        for (i.alpha in 1:n.alpha) {
          if (is.na(K.est[i.alpha]) & max(test.stat.K0$test.stat) <= threshold[i.alpha]) {
            K.est[i.alpha] <- K0
          }
          if (is.na(K.est.boot[i.alpha]) & max(test.stat.K0$test.stat.boot) <= threshold[i.alpha]) {
            K.est.boot[i.alpha] <- K0
          }
          if (!any(is.na(K.est)) & !any(is.na(K.est.boot))) {break}
        }
      }
    } else {
      for (K0 in K.min:K.max) {
        test.stat.K0 <- GoFTestBoot(A, K0 = K0, n.dim = n.dim, n.b = n.b, use.boot = use.boot)
        for (i.alpha in 1:n.alpha) {
          if (is.na(K.est[i.alpha]) & max(test.stat.K0$test.stat) <= threshold[i.alpha]) {
            K.est[i.alpha] <- K0
          }
          if (!any(is.na(K.est))) {break}
        }
      }    
    }
    return(list(K.est = K.est, K.est.boot = K.est.boot)) 
  }
  
  Generate.B <- function(K, low.limit = 0.1, upp.limit = 0.9) {
    B <- matrix(low.limit+ (upp.limit-low.limit)*runif(K^2), ncol = K)
    tB <- t(B)
    B[lower.tri(B)] <- tB[lower.tri(tB)]
    return(B)
  }
  
  Generate.clust.size <- function(n, K, prob = NULL){
    if (is.null(prob)) {
      prob <- rep(1/K, K)
    }
    clust.size <- rmultinom(1, n, prob)
    return(clust.size)
  }
  
  Generate.theta <- function(n, K, clust.size){
    theta <- matrix(0, n, K)
    for (k in 1:K){
      if (k==1){id1 <-1} else {id1 <- sum(clust.size[1:(k-1)])+1}
      id2 <- sum(clust.size[1:k])
      theta[id1:id2, k] <- 1
    } 
    return(theta)
  }
  
  Generate.psi <- function(n, K, clust.size, DCBM, low.limit = 0.2, upp.limit = 1){
    if (DCBM) {
      psi <- low.limit + (upp.limit-low.limit)* runif(n)
      for (k in 1:K){
        if (k==1){id1 <-1} else {id1 <- sum(clust.size[1:(k-1)])+1}
        id2 <- sum(clust.size[1:k])
        psi[id1:id2] <- psi[id1:id2] / max(psi[id1:id2])
      }
    } else {
      psi <- rep(1, n)
    }
    return(psi)
  }
  
  Generate.P <- function(n, K, model = 'sbm', dcbm.psi.bound = c(0, 1),
                         dirichlet.alpha = rep(1, K), threshold = 0.1,
                         low.limit = 0, upp.limit = 0.5) {
    if (model == 'sbm') {
      clust.size <- Generate.clust.size(n, K)
      Theta <- Generate.theta(n, K, clust.size)
      B <- Generate.B(K, low.limit = low.limit, upp.limit = upp.limit)
      while(svd(B)$d[K] < threshold) {
        B <- Generate.B(K, low.limit = low.limit, upp.limit = upp.limit)
      }
    }
    if (model == 'sbm.finer') {
      clust.size <- Generate.clust.size(n, K + 1)
      Theta <- Generate.theta(n, K + 1, clust.size)
      B <- Generate.B(K + 1, low.limit = low.limit, upp.limit = upp.limit)
      while(svd(B)$d[K + 1] < threshold) {
        B <- Generate.B(K + 1, low.limit = low.limit, upp.limit = upp.limit)
      }      
    }
    if (model == 'dcbm') {
      clust.size <- Generate.clust.size(n, K)
      psi <- Generate.psi(n, K, clust.size, DCBM = TRUE,
                          low.limit = dcbm.psi.bound[1], upp.limit = dcbm.psi.bound[2])
      Theta <- diag(psi) %*% Generate.theta(n, K, clust.size)
      B <- Generate.B(K, low.limit = low.limit, upp.limit = upp.limit)
      while(svd(B)$d[K] < threshold) {
        B <- Generate.B(K, low.limit = low.limit, upp.limit = upp.limit)
      }
    }
    if (model == 'mixed') {
      Theta <- rdirichlet(n, dirichlet.alpha)
      B <- Generate.B(K, low.limit = low.limit, upp.limit = upp.limit)
      while(svd(B)$d[K] < threshold) {
        B <- Generate.B(K, low.limit = low.limit, upp.limit = upp.limit)
      }
    }
    return(Theta %*% B %*% t(Theta))
  }
  
  Generate.P.diag <- function(n, K, model = 'sbm', dcbm.psi.bound = c(0, 1),
                              dirichlet.alpha = rep(1, K),
                              off.val = 0.1, diag.lift = 0.2) {
    if (model == 'sbm') {
      clust.size <- Generate.clust.size(n, K)
      Theta <- Generate.theta(n, K, clust.size)
      B <- matrix(off.val, ncol = K, nrow = K) + diag(rep(diag.lift, K))
    }
    if (model == 'sbm.finer') {
      clust.size <- Generate.clust.size(n, K + 1)
      Theta <- Generate.theta(n, K + 1, clust.size)
      B <- matrix(off.val, ncol = K + 1, nrow = K + 1) + diag(rep(diag.lift, K + 1)) 
    }
    if (model == 'dcbm') {
      clust.size <- Generate.clust.size(n, K)
      psi <- Generate.psi(n, K, clust.size, DCBM = TRUE,
                          low.limit = dcbm.psi.bound[1], upp.limit = dcbm.psi.bound[2])
      Theta <- diag(psi) %*% Generate.theta(n, K, clust.size)
      B <- matrix(off.val, ncol = K, nrow = K) + diag(rep(diag.lift, K))
    }
    if (model == 'mixed') {
      Theta <- rdirichlet(n, dirichlet.alpha)
      B <- matrix(off.val, ncol = K, nrow = K) + diag(rep(diag.lift, K))
    }
    return(Theta %*% B %*% t(Theta))
  }}
  




#n <- 500
#q <- c(1/3,1/3,1/3)

#a <- 0.3 #prob dentro da comunidade
#b <- 0.2 #prob fora das comunidades

#P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
#print(P)

#vetor de comunidades verdadeiro
#z0 <- sample.z(n, q)

#gerar matriz de adjacência

#X <- sample.g(z0,P)

#res <- SeqTestSBM(X)
#print(res)