library(mvtnorm)

interval_data_one = function(igama =  c(0.1,0.2)){
  K = 4
  p = 3
  nGrupos = rep(50,K)
  
  mu1 = c(-0.5,-0.5,0); sigma1 = diag(c(0.04,0.04,1),3) # Grupo 1
  mu2 = c(0.5,-0.5,0); sigma2 =  diag(c(0.04,0.04,1),3) # Grupo 2
  mu3 = c(0,0.5,-0.5); sigma3 =  diag(c(1,0.04,0.04),3) # Grupo 3
  mu4 = c(0,0.5,0.5); sigma4 =  diag(c(1,0.04,0.04),3)  # Grupo 4
  
  MU = list(mu1,mu2,mu3,mu4)
  SIGMA = array(NA,c(p,p,K))
  SIGMA[,,1] = sigma1; SIGMA[,,2] = sigma2; SIGMA[,,3] = sigma3; SIGMA[,,4] = sigma4
  
  number_list = list()
  group = NULL
  for(i in 1:K){
    number_list[[i]] = rmvnorm(nGrupos[i], mean = MU[[i]],sigma = SIGMA[,,i])
    group = c(group,rep(i,nGrupos[i]))
  }
  
  end_matrix = number_list[[1]]
  for(j in 2 : K)
    end_matrix = rbind( end_matrix, number_list[[j]])
  
  ## Creating an interval dataset 
  
  N = sum(nGrupos)
  interval_dset = matrix(NA,N,2*p)
  names_dset = NULL
  for(x in 1:p){
    l = runif(N,min = igama[1],max = igama[2])
    interval_dset[,2*x - 1] = end_matrix[,x] - l/2
    interval_dset[,2*x ] = end_matrix[,x] + l/2
    names_dset = c(names_dset, c(paste0("X",x,"_lower"),paste0("X",x,"_upper")))
  }
  interval_dset = cbind(interval_dset,group)
  colnames(interval_dset) = c(names_dset,"group")
  return(as.data.frame(interval_dset))
}

interval_data_two = function(igama = c(0.1,0.2)){
  K = 3
  p = 2
  nK = rep(30,3)
  
  mu1 = c(0.0, 0.0); sigma1 = diag(c(0.04,0.04),2)
  mu2 = c(0.1, 0.7); sigma2 = diag(c(0.04,0.04),2)
  mu3 = c(0.9, 0.1); sigma3 = diag(c(0.04,0.04),2)
  MU = list(mu1,mu2,mu3)
  SIGMA = array(NA,c(2,2,3))
  SIGMA[,,1] = sigma1; SIGMA[,,2] = sigma2; SIGMA[,,3] = sigma3
  
  N = sum(nK)
  number_list = list()
  group = NULL
  for(i in 1:K){
    
    number_list[[i]] = rmvnorm(nK[i], mean = MU[[i]],sigma = SIGMA[,,i])
    group = c(group,rep(i,nK[i]))
  }
  
  end_matrix = number_list[[1]]
  for(j in 2 : K)
    end_matrix = rbind( end_matrix, number_list[[j]])
  
  x3 = 2 * end_matrix[,1] - 1.5 * end_matrix[,2] + rnorm( N,0,1)
  
  end_matrix = cbind(end_matrix,x3)
  
  ## Creating an interval dataset 
  
  
  interval_dset = matrix(NA,N,2*(p+1) )
  names_dset = NULL
  for(x in 1:(p+1)){
    l = runif(N,min = igama[1],max = igama[2])
    interval_dset[,2*x - 1] = end_matrix[,x] - l/2
    interval_dset[,2*x ] = end_matrix[,x] + l/2
    names_dset = c(names_dset, c(paste0("X",x,"_lower"),paste0("X",x,"_upper")))
  }
  interval_dset = cbind(interval_dset,group)
  colnames(interval_dset) = c(names_dset,"group")
  return(as.data.frame(interval_dset))
}


