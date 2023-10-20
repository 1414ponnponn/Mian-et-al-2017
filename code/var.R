library(tidyverse)
library(pracma)

main_VAR <- function(data, index_y, index_x, p = 5, irhor = 11){
  result <- c()
  result$index_y <- index_y
  result$index_x <- index_x
  
  result$Y <- data[index_y] %>% 
    as.matrix()
  result$X <- data[index_x] %>% 
    as.matrix()
  result$obs <- nrow(result$Y)
  result$Ny <- ncol(result$Y)
  result$Nx <- ncol(result$X)
  
  bet <- c()
  res <- c()
  for(i in 1:result$Ny){
    y <- result$Y[,i]
    result_i <- lm(y ~ result$X)
    bet <- cbind(bet, result_i$coefficients)
    res <- cbind(res, result_i$residuals)
  }

  result$bet <- rbind(bet[2:nrow(bet),], bet[1,])
  result$res <- res
  result$Sigma <- t(res) %*% res / (result$obs - result$Nx)
  result$bet_vc <- matrix(result$bet, ncol = 1)
  
  result$ii <- histc(data[,1], unique(data[,1]))$cnt
  result$Ti <- histc(data[,1], unique(data[,1]))$bin
  
  result$Year <- data[,2]
  
  result$p <- p
  result$irhor <- irhor
  result$h <- 1:(irhor - 1)
  result$N <- max(result$ii)
  result$Perm <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
  result$normalize <- 1
  
  
  return(result)
}



correct_bias <- function(var_model){
  D <- zeros(result$obs, result$N)
  t <- 1
  for(i in 1:result$N){
    D[t:(t + result$Ti[i]-1), i] <- ones(result$Tie[i], 1)
    t <- t + result$Ti[i]
  }
  
  
}

irfvar <- function(var_model){
  Ny <- var_model$Ny
  p <- var_model$p
  Perm <- var_model$Perm
  Sigma <- var_model$Sigma
  bet <- var_model$bet
  nirf <- var_model$irhor
  
  Sig_chol <- t(t(Perm) %*% chol(Perm %*% Sigma %*% t(Perm)) %*% t(Perm))
  
  B_comp <- rbind(t(bet[1:(p*Ny),]), cbind(diag(nrow = Ny * (p-1)), zeros(Ny*(p-1), Ny)))
  
  cc_dy <- diag(Ny*p)
  msel <- rbind(diag(Ny), zeros((p-1) * Ny, Ny))
  
  if(var_model$normalize == 1){
    diagonal <- diag(diag(Sig_chol))
    Sig_chol <- Sig_chol %*% diag(diag(diagonal^-1))
  }
  
  virf <- array(0, dim = c(nirf, Ny, Ny))
  virf[1, ,] <- t(msel) %*% cc_dy %*% msel %*% Sig_chol
  
  for(tt in 2:nirf){
    cc_dy <- B_comp %*% cc_dy
    virf[tt, ,] <- t(msel) %*% cc_dy %*% msel %*% Sig_chol
  }

  var_model$Sig_chol <- Sig_chol
  var_model$cholIRF <- virf

  return(var_model)
}

chol_VAR <- function(var_model){
  h <- var_model$h
  
  varbc <- c()
  varbc$Ny <- var_model$Ny
  varbc$p <- var_model$p
  varbc$Perm <- var_model$Perm
  varbc$irhor <- var_model$irhor
  varbc$normalize <- var_model$normalize
  varbc$bet <- var_model$bet
  varbc$Sigma <- var_model$Sigma ##
  var_chol_bc <- irfvar(varbc)
}

plot_irf <- function(vec_irf){
  len <- length(vec_irf)
  g <- data.frame(t = 0:(len-1), irf = vec_irf) %>% 
    ggplot(aes(x = t, y = irf)) +
    geom_point(shape = 8, col = "red") +
    geom_line(col = "red")
}
