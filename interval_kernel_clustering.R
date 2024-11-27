# All codes reference the cluster methods studied in my scientific initiation project

# Fun??es usadas nos m?todos de agrupamento
library(rdist)

sigma2est = function(X) {
  X <- as.matrix(X)
  h = 1:ncol(X)
  xmin <- as.matrix(X[,which(h%%2==1)])
  xmax <- as.matrix(X[,which(h%%2==0)])
  d_min = as.vector(proxy::dist(xmin))**2
  d_max = as.vector(proxy::dist(xmax))**2
  dist = d_min + d_max
  mean(quantile(dist[dist != 0], probs = c(.9, .1)))
}



# Helper functions to compute the kernel matrices

itrval_gaussian_kernel <- function(A_lower, A_upper, B_lower, B_upper, sig2){
  D2 <-   (A_lower - B_lower)**2 + (A_upper - B_upper)**2
  const <- 1/(2*sig2)
  return( exp(-const * D2) )
}

single_kernel_matrix <- function (Xi, Xs, Yi, Ys, sig2){
  n1 = length(Xi); n2 = length(Yi)
  Xi <- rep(Xi, rep.int(length(Yi), length(Xi)))
  Xs <- rep(Xs, rep.int(length(Ys), length(Xs)))
  Yi <- rep(Yi, times = ceiling(length(Xi)/length(Yi)))
  Ys <- rep(Ys, times = ceiling(length(Xs)/length(Ys)))
  MK <- itrval_gaussian_kernel(Xi, Xs, Yi, Ys, sig2)
  MK <- matrix(MK,n1,n2,byrow = T)
  return(MK)
}

kernel_matrix <- function(A, B, Ga, Gb, sig2, n, p, K){
  KMs <- array(NA,c(n,K,p))
  for(j in 1:p)
    KMs[,,j] <- single_kernel_matrix(A[,j], B[,j], Ga[,j], Gb[,j],sig2)
  return(KMs)
}

matrix_kernel_conventional <- function(A, B, sig2){
  D2a = rdist::cdist(X = A, Y = A, metric = 'euclidean', p = 2L) ** 2
  D2b = rdist::cdist(X = B, Y = B, metric = 'euclidean', p = 2L) ** 2
  D2 = D2a + D2b
  const = -1/(2*sig2)
  return( exp( const * D2) )
}



# Helper functions to compute the distances

dist_in_kernel_space_per_variable <- function(A, B, Part, sig2, n, p, K, d2 = NULL){
  
  if(is.null(d2)){
    d2 <- array(Inf,c(n,K,p))
  }
  
  for(k in 1:K){
    P <- which( Part == k)
    card <- length(P)
    if(card == 0){
      next
    }
    else{
      for(j in 1:p){
        Aj = matrix(A[,j], ncol = 1)
        Bj = matrix(B[,j], ncol = 1)
        KM = matrix_kernel_conventional(A = Aj, B = Bj, sig2 = sig2)
        KMk = KM[P,P]
        second_sum = sum(KMk)/(card**2)
        if(card == 1){
          first_sum = (2/card) * KM[,P]
        }
        else{
          first_sum = (2/card) * apply(X = KM[,P], MARGIN = 1, FUN = sum)
        }
        d2[,k,j] = 1.0 - first_sum + second_sum
      }
    }
  }
  return(d2)
}

d2_weighted <- function(d2j, var_weights){
  
  n = dim(d2j)[1]; K= dim(d2j)[2]; p = dim(d2j)[3]
  D = matrix(Inf, nrow = n, ncol = K)
  
  if( is.vector(var_weights) ){
    var_weights = matrix(rep(var_weights, K), nrow = K, byrow = TRUE)
  }
  
  for(k in 1:K){
    D[,k] = d2j[,k,]  %*%  var_weights[k,]
  }
  return(D)
}




# Helper functions for computing the weights of the variables

get_weights_GP <- function(d2j, Part, K, p, var_weights ){
  Dkj <- matrix(NA, K, p)
  for (k in 1:K) {
    P <- which(Part == k)
    nk = length(P)
    if( nk == 0){
      Dkj[k,] = 0.0
    }
    else{
      if(nk==1){
        Dkj[k,] = d2j[P,k,]
      }
      else{
        Dkj[k,] = apply(X = d2j[P,k,], MARGIN = 2, FUN = sum)
      }
    }
  }
  Dj = apply(X = Dkj, MARGIN = 2, FUN = sum)
  if(all(Dj > 1e-20)){
    var_weights <- prod(Dj ** (1/p))/Dj
    return(var_weights)
  }
  else{
    idx_lower = which(Dj <= 1e-20)
    idx_greater = which(Dj > 1e-20)
    n_lower = length(idx_lower)
    
    if(n_lower == p){
      return(var_weights)
    }
    else{
      const = (1/prod(var_weights[idx_lower]))** (1/(p-n_lower))
      Dj_greater = Dj[idx_greater]
      var_weights[idx_greater] = const * (prod(Dj_greater ** (1/p))/Dj_greater)
      return(var_weights)
    }
  }
}

get_weights_LP <- function(d2j, Part, K, p, var_weights ){
  for (k in 1:K) {
    P <- which(Part == k)
    nk = length(P)
    if( nk == 0){
      next
    }
    else{
      if(nk==1){
        Dj = d2j[P,k,]
      }
      else{
        Dj = apply(X = d2j[P,k,], MARGIN = 2, FUN = sum)
      }
    }

    if(all(Dj > 1e-20)){
      var_weights[k,] <- prod(Dj ** (1/p))/Dj
    }
    else{
      idx_lower = which(Dj <= 1e-20)
      idx_greater = which(Dj > 1e-20)
      n_lower = length(idx_lower)
      
      if(n_lower == p){
        next
      }
      else{
        const = (1/prod(var_weights[k, idx_lower]))** (1/(p-n_lower))
        Dj_greater = Dj[idx_greater]
        var_weights[k, idx_greater] = const * (prod(Dj_greater ** (1/p))/Dj_greater)
      }
    }
  }
  return(var_weights)
}

get_weights_GS <- function(d2j, theta, Part, K, p, var_weights){
  exponent = 1/(theta-1)
  
  Dkj <- matrix(NA, K, p)
  for (k in 1:K) {
    P <- which(Part == k)
    nk = length(P)
    if( nk == 0){
      Dkj[k,] = 0.0
    }
    else{
      if(nk==1){
        Dkj[k,] = d2j[P,k,]
      }
      else{
        Dkj[k,] = apply(X = d2j[P,k,], MARGIN = 2, FUN = sum)
      }
    }
  }
  Dj = apply(X = Dkj, MARGIN = 2, FUN = sum)
  inv_dj = (1/Dj) ** exponent
  idx_inf = which(is.infinite(inv_dj))
  n_inf = length(idx_inf)
  
  if(n_inf == 0){
    return(inv_dj/sum(inv_dj))
  }
  else{
    if(n_inf == p){
      return(var_weights)
    }
    else{
      idx_non = which(!is.infinite(inv_dj))
      inv_dj_new = inv_dj[idx_non]
      sum_previous = sum( var_weights[idx_inf] )
      const = 1 - sum_previous
      var_weights[idx_non] = const * (inv_dj_new/sum(inv_dj_new))
    }
  }
  return(var_weights)
}

get_weights_LS <- function(d2j, theta, Part, K, p, var_weights){
  exponent = 1/(theta-1)
  for (k in 1:K) {
    P <- which(Part == k)
    nk = length(P)
    if( nk == 0){
      next
    }
    else{
      if(nk==1){
        Dj = d2j[P,k,]
      }
      else{
        Dj = apply(X = d2j[P,k,], MARGIN = 2, FUN = sum)
      }
    }
    
    inv_dj = (1/Dj) ** exponent
    idx_inf = which(is.infinite(inv_dj))
    n_inf = length(idx_inf)
    
    if(n_inf == 0){
      var_weights[k,] = inv_dj/sum(inv_dj)
    }
    else{
      if(n_inf == p){
        next
      }
      else{
        idx_non = which(!is.infinite(inv_dj))
        inv_dj_new = inv_dj[idx_non]
        sum_previous = sum( var_weights[k, idx_inf] )
        const = 1 - sum_previous
        var_weights[k, idx_non] = const * (inv_dj_new/sum(inv_dj_new))
      }
    }
  }
  return(var_weights)
}


# Helper functions for computing the prototypes 

initial_prototypes = function(X, K){
  n = nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  xx <- unique(X); nn <- nrow(xx)
  is <- sample.int(nn,K)
  return( list(lower_prototypes = xx[is, idl], upper_prototypes = xx[is, idu]) )
}

prototypes_kernel <- function(A, B, KM, Part, K, Ga, Gb){
  ld = 1e-10
  for(k in 1:K){
    P = which(Part == k)
    nk = length(P)
    if(nk == 0){
      next
    }
    
    if(nk == 1){
      KMj = KM[P,k,]
      idx_upper = which(KMj > ld)
      ga = (KMj * A[P,])/KMj
      gb = (KMj * B[P,])/KMj
      Ga[k,idx_upper] = ga[idx_upper]
      Gb[k,idx_upper] = gb[idx_upper]
    }
    else{ 
      KMj = apply(X = KM[P,k,], MARGIN = 2, FUN = sum)
      idx_upper = which(KMj > ld)
      numerator_lower = apply(X = KM[P,k,] * A[P,], MARGIN = 2, FUN = sum)
      numerator_upper = apply(X = KM[P,k,] * B[P,], MARGIN = 2, FUN = sum)
      ga = numerator_lower/KMj
      gb = numerator_upper/KMj
      Ga[k,idx_upper] = ga[idx_upper]
      Gb[k,idx_upper] = gb[idx_upper]
      
    }
  }
  return(list(lower_prototypes = Ga, upper_prototypes = Gb))
}



# Helper function for computing the objective function

get_J_general <- function(D, P, K){
  J <- 0
  for(k in 1:K){
    if(k %in% P){
      J <- J + sum( D[P == k,k] )
    }
  }
  return(J)
}



##########################################################################################################
###                                         Main Functions                                             ###
##########################################################################################################



IVKKM_O <- function(X, K, sig2, maxiter = 100){
  
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  G = initial_prototypes(X = X, K = K)
  Ga = G$lower_prototypes; Gb = G$upper_prototypes
  
  KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
  d2_per_variable = 2-2*KM
  d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
  
  P <- apply(d2,1, which.min)
  J <-  get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # updating the prototypes
    G <- prototypes_kernel(A = A, B = B, KM = KM, Part = P, K = K, Ga = Ga, Gb = Gb)
    Ga = G$lower_prototypes; Gb = G$upper_prototypes
    
    # computing the kernel matrix fo each variable
    KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
    
    # computing the distances
    d2_per_variable = 2-2*KM
    d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = P, K = K)
    }
    
    J <- c(J,Jn)
  }
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J))
    warning("The objective function did not decrease")
  
  return( list(partition = P,lower_centroids = Ga, upper_centroids = Gb, iterations = it, criterion = J) )
  
}

IVKKM_O_GP <- function(X, K, sig2, maxiter = 100){
  
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights = rep(1, p)
  
  G = initial_prototypes(X = X, K = K)
  Ga = G$lower_prototypes; Gb = G$upper_prototypes
  
  KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
  d2_per_variable = 2-2*KM
  d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
  
  P <- apply(d2,1, which.min)
  J <-  get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # updating the prototypes
    G <- prototypes_kernel(A = A, B = B, KM = KM, Part = P, K = K, Ga = Ga, Gb = Gb)
    Ga = G$lower_prototypes; Gb = G$upper_prototypes
    
    # computing the kernel matrix fo each variable
    KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
    
    # computing the distances per variable
    d2_per_variable = 2-2*KM
    
    # computing the weights of the variables
    var_weights = get_weights_GP(d2j = d2_per_variable, Part = P, K = K, p = p, var_weights = var_weights)
    
    # computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = var_weights)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = P, K = K)
    }
    
    J <- c(J,Jn)
  }
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J))
    warning("The objective function did not decrease")
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  names(var_weights) <- name.var
  
  return( list(partition = P,lower_centroids = Ga, upper_centroids = Gb, variable_weights = var_weights, iterations = it, criterion = J) )
}

IVKKM_O_LP <- function(X, K, sig2, maxiter = 100){
  
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights = matrix(data = 1, nrow = K, ncol = p)
  
  G = initial_prototypes(X = X, K = K)
  Ga = G$lower_prototypes; Gb = G$upper_prototypes
  
  KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
  d2_per_variable = 2-2*KM
  d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
  
  P <- apply(d2,1, which.min)
  J <-  get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # updating the prototypes
    G <- prototypes_kernel(A = A, B = B, KM = KM, Part = P, K = K, Ga = Ga, Gb = Gb)
    Ga = G$lower_prototypes; Gb = G$upper_prototypes
    
    # computing the kernel matrix fo each variable
    KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
    
    # computing the distances per variable
    d2_per_variable = 2-2*KM
    
    # computing the weights of the variables
    var_weights = get_weights_LP(d2j = d2_per_variable, Part = P, K = K, p = p, var_weights = var_weights)
    
    # computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = var_weights)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = P, K = K)
    }
    
    J <- c(J,Jn)
  }
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J))
    warning("The objective function did not decrease")
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  colnames(var_weights) <- name.var
  
  return( list(partition = P,lower_centroids = Ga, upper_centroids = Gb, variable_weights = var_weights, iterations = it, criterion = J) )
}

IVKKM_O_GS <- function(X, K, theta, sig2, maxiter = 100){
  
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights = rep(1/p, p)
  
  G = initial_prototypes(X = X, K = K)
  Ga = G$lower_prototypes; Gb = G$upper_prototypes
  
  KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
  d2_per_variable = 2-2*KM
  d2 = d2_weighted(d2j = d2_per_variable, var_weights = (var_weights**theta) )
  
  P <- apply(d2,1, which.min)
  J <-  get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # updating the prototypes
    G <- prototypes_kernel(A = A, B = B, KM = KM, Part = P, K = K, Ga = Ga, Gb = Gb)
    Ga = G$lower_prototypes; Gb = G$upper_prototypes
    
    # computing the kernel matrix fo each variable
    KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
    
    # computing the distances per variable
    d2_per_variable = 2-2*KM
    
    # computing the weights of the variables
    var_weights = get_weights_GS(d2j = d2_per_variable, theta = theta, Part = P, K = K, p = p, var_weights = var_weights)
    
    # computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = var_weights**theta)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = P, K = K)
    }
    
    J <- c(J,Jn)
  }
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J))
    warning("The objective function did not decrease")
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  names(var_weights) <- name.var
  
  return( list(partition = P,lower_centroids = Ga, upper_centroids = Gb, variable_weights = var_weights, iterations = it, criterion = J) )
}

IVKKM_O_LS <- function(X, K, theta, sig2, maxiter = 100){
  
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights = matrix(data = 1/p, nrow = K, ncol = p)
  
  G = initial_prototypes(X = X, K = K)
  Ga = G$lower_prototypes; Gb = G$upper_prototypes
  
  KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
  d2_per_variable = 2-2*KM
  d2 = d2_weighted(d2j = d2_per_variable, var_weights = (var_weights**theta) )
  
  P <- apply(d2,1, which.min)
  J <-  get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # updating the prototypes
    G <- prototypes_kernel(A = A, B = B, KM = KM, Part = P, K = K, Ga = Ga, Gb = Gb)
    Ga = G$lower_prototypes; Gb = G$upper_prototypes
    
    # computing the kernel matrix fo each variable
    KM = kernel_matrix(A = A, B = B, Ga = Ga, Gb = Gb, sig2 = sig2, n = n, p = p, K = K)
    
    # computing the distances per variable
    d2_per_variable = 2-2*KM
    
    # computing the weights of the variables
    var_weights = get_weights_LS(d2j = d2_per_variable, theta = theta, Part = P, K = K, p = p, var_weights = var_weights)
    
    # computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = var_weights**theta)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = P, K = K)
    }
    
    J <- c(J,Jn)
  }
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J))
    warning("The objective function did not decrease")
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  colnames(var_weights) <- name.var
  
  return( list(partition = P,lower_centroids = Ga, upper_centroids = Gb, variable_weights = var_weights, iterations = it, criterion = J) )
}



IVKKM_K <- function(X, K, sig2, maxiter = 100){
  
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  # Initial partition
  P <- sample(K,n,replace = TRUE)
  
  # computing the distances per variable
  d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = NULL)
  
  # Computing the final distance
  d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
  
  # Updating the partition
  P <- apply(d2, 1, which.min)
  J = get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    
    it <- it + 1
    # computing the distances per variable
    d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = d2_per_variable)
    
    # Computing the final distance
    d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = P, K = K)
    }
    J <- c(J,Jn)
  }
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J) )
    warning("The objective function did not decrease")
  
  return( list(partition = P, iterations = it, criterion = J) )
}

IVKKM_K_GP <- function(X, K, sig2, maxiter = 100){
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights <- rep(1,p)
  
  # Initial partition
  P <- sample(K,n,replace = TRUE)
  
  # Computing the distance per variable
  d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = NULL)
  
  # Computing the final distance
  d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
  
  # Updating the partition
  P <- apply(d2, 1, which.min)
  
  
  J = get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # Computing the distance per variable
    d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = d2_per_variable)
    
    # Computing the weights of the variables
    var_weights = get_weights_GP(d2j = d2_per_variable, Part = P, K = K, p = p, var_weights = var_weights)
    
    # Computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = var_weights)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = Pn, K = K)
    }
    J <- c(J,Jn)
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  names(var_weights) <- name.var
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J) )
    warning("The objective function did not decrease")
  
  return( list(partition = P, variable_weights = var_weights, iterations = it, criterion = J) )
  
}

IVKKM_K_GS <- function(X, K, theta, sig2, maxiter = 100){
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights <- rep(1/p, p)
  
  # Initial partition
  P <- sample(K,n,replace = TRUE)
  
  # Computing the distance per variable
  d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = NULL)
  

  # Computing the final distance
  d2 = d2_weighted(d2j = d2_per_variable, var_weights = (var_weights**theta) )
  
  # Updating the partition
  P <- apply(d2, 1, which.min)
  
  J = get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # Computing the distance per variable
    d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = d2_per_variable)
    
    # Computing the weights of the variables
    var_weights = get_weights_GS(d2j = d2_per_variable, theta = theta, Part = P, K = K, p = p, var_weights = var_weights)
    
    # Computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = (var_weights**theta) )
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = Pn, K = K)
    }
    J <- c(J,Jn)
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  names(var_weights) <- name.var
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J) )
    warning("The objective function did not decrease")
  
  return( list(partition = P, variable_weights = var_weights, iterations = it, criterion = J) )
  
}

IVKKM_K_LP <- function(X, K, sig2, maxiter = 100){
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights <- matrix(data = 1, nrow = K, ncol = p)
  
  # Initial partition
  P <- sample(K,n,replace = TRUE)
  
  # Computing the distance per variable
  d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = NULL)
  
  # Computing the final distance
  d2 = apply(X = d2_per_variable, MARGIN = c(1,2), FUN = sum)
  
  # Updating the partition
  P <- apply(d2, 1, which.min)
  
  J = get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # Computing the distance per variable
    d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = d2_per_variable)
    
    # Computing the weights of the variables
    var_weights = get_weights_LP(d2j = d2_per_variable, Part = P, K = K, p = p, var_weights = var_weights)
    
    # Computing the final distance
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = var_weights)
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = Pn, K = K)
    }
    J <- c(J,Jn)
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  colnames(var_weights) <- name.var
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J) )
    warning("The objective function did not decrease")
  
  return( list(partition = P, variable_weights = var_weights, iterations = it, criterion = J) )
  
}

IVKKM_K_LS <- function(X, K, theta, sig2, maxiter = 100){
  X <- as.matrix(X)
  n <- nrow(X)
  nc = ncol(X); idl = seq(1,nc, by = 2); idu = seq(2,nc, by = 2)
  A = as.matrix(X[, idl])
  B = as.matrix(X[, idu])
  p = ncol(A)
  
  var_weights <- matrix(data = 1/p, nrow = K, ncol = p)
  
  # Initial partition
  P <- sample(K,n,replace = TRUE)
  
  # Computing the distance per variable
  d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = NULL)
  
  # Computing the final distance
  d2 = d2_weighted(d2j = d2_per_variable, var_weights = (var_weights**theta) )
  
  # Updating the partition
  P <- apply(d2, 1, which.min)
  
  J = get_J_general(D = d2, P = P, K = K)
  
  it <- 0
  repeat{
    it <- it + 1
    
    # Computing the distance per variable
    d2_per_variable = dist_in_kernel_space_per_variable(A = A, B = B, Part = P, sig2 = sig2, n = n, p = p, K = K, d2 = d2_per_variable)
    
    # Computing the weights of the variables
    var_weights = get_weights_LS(d2j = d2_per_variable, theta = theta, Part = P, K = K, p = p, var_weights = var_weights)
    
    # Computing the weights of the variables
    d2 = d2_weighted(d2j = d2_per_variable, var_weights = (var_weights**theta) )
    
    # Updating the partition
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || it >= maxiter) break
    else{
      P <- Pn
      Jn = get_J_general(D = d2, P = Pn, K = K)
    }
    J <- c(J,Jn)
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[seq(1,nc, by = 2)]
  colnames(var_weights) <- name.var
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( !all(Jsort == J) )
    warning("The objective function did not decrease")
  
  return( list(partition = P, variable_weights = var_weights, iterations = it, criterion = J) )
  
}




