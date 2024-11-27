source("path_to_the codes.R")
library(ggplot2)


# Function to standardize the data according to the midpoints of the intervals
standarlization <- function(X){
  X <- as.matrix(X)
  n = nrow(X); d = ncol(X)
  idl = seq(1,d, by = 2); idu = seq(2,d, by = 2)
  A = as.matrix(X[, idl]); B = as.matrix(X[, idu])
  p = ncol(A)
  
  Mi = (A + B)/2
  M = apply(X = Mi, MARGIN = 2, FUN = mean)
  
  Si = t( t(Mi) - M) ** 2
  S =  apply(X = Si, MARGIN = 2, FUN = mean)
  
  Dl = t( t(A) - M)
  Du = t( t(B) - M)
  
  A_new = t(t(Dl)/sqrt(S))
  B_new = t(t(Du)/sqrt(S))
  
  X_new = matrix(NA, nrow = n, ncol = d)
  X_new[,idl] = A_new
  X_new[,idu] = B_new
  colnames(X_new) = colnames(X)
  return(X_new)
}



#####################################################
###       Application to the Fungi Dataset        ###
#####################################################


dataset  = read.table("path_to_fungi_dataset.txt", header = TRUE)

y = dataset[, 11] # selecting the labels of the dataset
X = standarlization(dataset[,-11]) # standardizing the dataset
sig2 = sigma2est(X) # Estimating the value of sigma^2

#----------------------------------#
# Applying the algorithm IVKKM-O   #
#----------------------------------#

set.seed(2)
rs = IVKKM_O(X = X, K = 3, sig2 = sig2)


# Obtaining the partition 
rs$partition

# Obtaining the lower centroids
rs$lower_centroids  # For the methods in kernel space there are no centroids

# Obtaining the upper centroids
rs$upper_centroids  # For the methods in kernel space there are no centroids

# Obtaining the weights of the variables
rs$criterion

data <- data.frame(iterations = 1:rs$iterations, J = rs$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-O", x = "Iterations", y = "Objective Function") +
  theme_minimal()                



#----------------------------------#
# Applying the algorithm IVKKM-O-LP   #
#----------------------------------#

set.seed(2)
rs_olp = IVKKM_O_LP(X = X, K = 3, sig2 = sig2)


# Obtaining the partition 
rs_olp$partition

# Obtaining the lower centroids
rs_olp$lower_centroids  # For the methods in kernel space there are no centroids

# Obtaining the upper centroids
rs_olp$upper_centroids  # For the methods in kernel space there are no centroids

# Obtaining the weights of the variables
rs_olp$variable_weights

# Obtaining the objective function 
rs_olp$criterion

data <- data.frame(iterations = 1:rs_olp$iterations, J = rs_olp$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-O-LP", x = "Iterations", y = "Objective Function") +
  theme_minimal()  


#----------------------------------#
# Applying the algorithm IVKKM-O-GP   #
#----------------------------------#

set.seed(2)
rs_ogp = IVKKM_O_GP(X = X, K = 3, sig2 = sig2)


# Obtaining the partition 
rs_ogp$partition

# Obtaining the lower centroids
rs_ogp$lower_centroids  # For the methods in kernel space there are no centroids

# Obtaining the upper centroids
rs_ogp$upper_centroids  # For the methods in kernel space there are no centroids

# Obtaining the weights of the variables
rs_ogp$variable_weights

# Obtaining the objective function 
rs_ogp$criterion

data <- data.frame(iterations = 1:rs_ogp$iterations, J = rs_ogp$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-O-GP", x = "Iterations", y = "Objective Function") +
  theme_minimal()  



#----------------------------------#
# Applying the algorithm IVKKM-O-LS   #
#----------------------------------#

set.seed(2)
rs_ols = IVKKM_O_LS(X = X, K = 3, theta = 1.5, sig2 = sig2)


# Obtaining the partition 
rs_ols$partition

# Obtaining the lower centroids
rs_ols$lower_centroids  # For the methods in kernel space there are no centroids

# Obtaining the upper centroids
rs_ols$upper_centroids  # For the methods in kernel space there are no centroids

# Obtaining the weights of the variables
rs_ols$variable_weights

# Obtaining the weights of the variables
rs_ols$criterion

data <- data.frame(iterations = 1:rs_ols$iterations, J = rs_ols$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-O-LS", x = "Iterations", y = "Objective Function") +
  theme_minimal()  



#----------------------------------#
# Applying the algorithm IVKKM-O-GS   #
#----------------------------------#

set.seed(2)
rs_ogs = IVKKM_O_GS(X = X, K = 3, theta = 1.5, sig2 = sig2)


# Obtaining the partition 
rs_ogs$partition

# Obtaining the lower centroids
rs_ogs$lower_centroids  # For the methods in kernel space there are no centroids

# Obtaining the upper centroids
rs_ogs$upper_centroids  # For the methods in kernel space there are no centroids

# Obtaining the weights of the variables
rs_ogs$variable_weights

# Obtaining the weights of the variables
rs_ogs$criterion

data <- data.frame(iterations = 1:rs_ogs$iterations, J = rs_ogs$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-O-GS", x = "Iterations", y = "Objective Function") +
  theme_minimal()  







#----------------------------------#
# Applying the algorithm IVKKM-K   #
#----------------------------------#

set.seed(2)
rs_k = IVKKM_K(X = X, K = 3, sig2 = sig2)


# Obtaining the partition 
rs_k$partition

# Obtaining the weights of the variables
rs_k$criterion

data <- data.frame(iterations = 1:rs_k$iterations, J = rs_k$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-K", x = "Iterations", y = "Objective Function") +
  theme_minimal()                



#----------------------------------#
# Applying the algorithm IVKKM-K-LP   #
#----------------------------------#

set.seed(2)
rs_klp = IVKKM_K_LP(X = X, K = 3, sig2 = sig2)


# Obtaining the partition 
rs_klp$partition

# Obtaining the weights of the variables
rs_klp$variable_weights

# Obtaining the objective function 
rs_klp$criterion

data <- data.frame(iterations = 1:rs_klp$iterations, J = rs_klp$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-K-LP", x = "Iterations", y = "Objective Function") +
  theme_minimal()  


#----------------------------------#
# Applying the algorithm IVKKM-K-GP   #
#----------------------------------#

set.seed(2)
rs_kgp = IVKKM_K_GP(X = X, K = 3, sig2 = sig2)


# Obtaining the partition 
rs_kgp$partition


# Obtaining the weights of the variables
rs_kgp$variable_weights

# Obtaining the objective function 
rs_kgp$criterion

data <- data.frame(iterations = 1:rs_kgp$iterations, J = rs_kgp$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-K-GP", x = "Iterations", y = "Objective Function") +
  theme_minimal()  



#----------------------------------#
# Applying the algorithm IVKKM-O-LS   #
#----------------------------------#

set.seed(2)
rs_kls = IVKKM_K_LS(X = X, K = 3, theta = 1.5, sig2 = sig2)


# Obtaining the partition 
rs_kls$partition

# Obtaining the weights of the variables
rs_kls$variable_weights

# Obtaining the weights of the variables
rs_kls$criterion

data <- data.frame(iterations = 1:rs_kls$iterations, J = rs_kls$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-K-LS", x = "Iterations", y = "Objective Function") +
  theme_minimal()  



#----------------------------------#
# Applying the algorithm IVKKM-K-GS   #
#----------------------------------#

set.seed(2)
rs_kgs = IVKKM_K_GS(X = X, K = 3, theta = 1.5, sig2 = sig2)


# Obtaining the partition 
rs_kgs$partition

# Obtaining the weights of the variables
rs_kgs$variable_weights

# Obtaining the weights of the variables
rs_kgs$criterion

data <- data.frame(iterations = 1:rs_kgs$iterations, J = rs_kgs$criterion)
ggplot(data, aes(x = iterations, y = J)) +
  geom_line() +    
  geom_point() +    
  labs(title = "IVKKM-K-GS", x = "Iterations", y = "Objective Function") +
  theme_minimal()  



