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
  result$Nx <- ncol(result$X) + 1
  
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
  result$bet_bc <- matrix(result$bet, ncol = 1)
  
  result$ii <- histc(data[,1], unique(data[,1]))$bin
  result$Ti <- histc(data[,1], unique(data[,1]))$cnt
  
  result$Year <- data[,2]
  
  result$p <- p
  result$irhor <- irhor
  result$h <- 1:(irhor - 1)
  result$N <- max(result$ii)
  result$Perm <- rbind(c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))
  result$normalize <- 1
  
  
  return(result)
}

doVAR_BiasCorrected <- function(VAR, nboot, clevel, tol, delta_update) {
  res <- detrend(VAR$res, tt = "constant")
  conv_test <- 1
  bet_bs <- VAR$bet0[1:(VAR$Ny * VAR$p), , drop = FALSE]
  
  while (conv_test > tol) {
    bet_bs_smp <- array(0, dim = c(VAR$Ny * VAR$p, VAR$Ny, nboot))
    
    for (jj in 1:nboot) {
      rr <- sample(c(-1, 1), size = VAR$obs, replace = TRUE)
      resb <- res * (rr %*% matrix(1, 1, VAR$Ny))
      t <- 1
      Yb <- matrix(0, nrow = VAR$obs, ncol = VAR$Ny)
      Xb <- matrix(0, nrow = VAR$obs, ncol = VAR$Ny * VAR$p)
      
      for (i in 1:VAR$N) {
        Xb[t, ] <- VAR$X[t, 1:(VAR$Ny * VAR$p)]
        for (j in t:(t + VAR$Ti[i] - 1)) {
          Yb[j, ] <- Xb[j, ] %*% bet_bs[1:(VAR$Ny*VAR$p),] + resb[j, , drop = FALSE]
          if (j < t + VAR$Ti[i] - 1) {
            Xb[j + 1, ] <- c(Yb[j, ], Xb[j, 1:((VAR$p - 1) * VAR$Ny)])
          }
        }
        t <- t + VAR$Ti[i]
      }
      bet_tmp <- qr.solve(cbind(Xb, VAR$D[, 1:(VAR$N - 1)], matrix(1, nrow(Xb), 1)),  Yb)
      bet_bs_smp[, , jj] <- bet_tmp[1:(VAR$Ny * VAR$p), , drop = FALSE]
    }
    
    betbs_mean <- apply(bet_bs_smp, c(1, 2), mean)
    diff <- VAR$bet0[1:(VAR$Ny * VAR$p), ] - betbs_mean[1:(VAR$Ny*VAR$p),]
    diffvec <- as.vector(diff)
    cat("Convergence test statistic (root mean squared difference) = ", sqrt(mean(diffvec^2)), "\n")
    conv_test <- sqrt(mean(diffvec^2))
    
    if (conv_test > tol) {
      bet_bs <- bet_bs + delta_update * diff
    }
  }
  
  VARbc <- list()
  VARbc$bet <- bet_bs
  
  return(VARbc)
}


correct_bias <- function(VAR){
  # Step 0: Create D (panel dummy variable matrix), demeaned Y, demeaned X
  D <- matrix(0, nrow = VAR$obs, ncol = VAR$N)
  t <- 1
  for (i in 1:VAR$N) {
    D[t:(t + VAR$Ti[i] - 1), i] <- rep(1, VAR$Ti[i])
    t <- t + VAR$Ti[i]
  }
  
  # Demean Y and X
  b_means <- solve(t(D) %*% D) %*% t(D) %*% cbind(VAR$Y, VAR$X[, 1:(VAR$Ny * VAR$p)])
  VAR$Ydm <- VAR$Y - D %*% b_means[, 1:VAR$Ny]
  VAR$Xdm <- VAR$X[, 1:(VAR$Ny * VAR$p)] - D %*% b_means[, (VAR$Ny + 1):ncol(b_means)]
  VAR$bettil <- qr.solve(VAR$Xdm, VAR$Ydm)
  VAR$D <- D
  
  # Step 1: Define inputs to bias correction function
  tol <- 0.0015
  nboot <- 2000
  delta_update <- 0.5
  VAR$bet0 <- VAR$bet  # Initial guess for bet_bc is OLS beta
  
  # Step 2: Run bias correction function
  VARbc <- doVAR_BiasCorrected(VAR, nboot, clevel, tol, delta_update)
  
  # Step 3: Extract bias and bias-corrected beta from VARbc
  VAR$bc_diff <- VAR$bet[1:(VAR$Ny * VAR$p), ] - VARbc$bet[1:(VAR$Ny * VAR$p), ]  # Bias
  VAR$bet_bc <- rbind(VARbc$bet, VAR$bet[(VAR$Ny * VAR$p + 1):nrow(VAR$bet), ])  # Bias-corrected beta
  
  # Adjust fixed effects
  uhat <- VAR$Y - VAR$X[, 1:(VAR$Ny * VAR$p)] %*% VAR$bet_bc[1:(VAR$Ny * VAR$p), ]
  D1 <- cbind(D[, 1:(VAR$N - 1)], rep(1, VAR$obs))
  alpha_bc <- solve(t(D1) %*% D1) %*% t(D1) %*% uhat
  VAR$bet_bc[(VAR$Ny * VAR$p + 1):nrow(VAR$bet_bc), ] <- alpha_bc
  VAR$res_bc <- VAR$Y - cbind(VAR$X, rep(1, VAR$obs)) %*% VAR$bet_bc
  VAR$Sigma_bc <- t(VAR$res_bc) %*% VAR$res_bc / (VAR$obs - VAR$Nx)  # Bias-corrected Sigma
  VAR$res <- VAR$res_bc  # Update residuals
  
  # Residuals need to be rescaled to account for the fact that estimated error terms have lower variance than population error terms
  VAR$res <- VAR$res * sqrt(VAR$obs / (VAR$obs - VAR$Nx))
  VAR$Sigma <- t(VAR$res) %*% VAR$res / (VAR$obs - VAR$Nx)
  
  return(VAR)
}

irfvar <- function(VAR){
  Ny <- VAR$Ny
  p <- VAR$p
  Perm <- VAR$Perm
  Sigma <- VAR$Sigma
  bet <- VAR$bet
  nirf <- VAR$irhor
  
  Sig_chol <- t(t(Perm) %*% chol(Perm %*% Sigma %*% t(Perm)) %*% t(Perm))
  
  B_comp <- rbind(t(bet[1:(p*Ny),]), cbind(diag(nrow = Ny * (p-1)), zeros(Ny*(p-1), Ny)))
  
  cc_dy <- diag(Ny*p)
  msel <- rbind(diag(Ny), zeros((p-1) * Ny, Ny))
  
  if(VAR$normalize == 1){
    diagonal <- diag(diag(Sig_chol))
    Sig_chol <- Sig_chol %*% diag(diag(diagonal^-1))
  }
  
  virf <- array(0, dim = c(nirf, Ny, Ny))
  virf[1, ,] <- t(msel) %*% cc_dy %*% msel %*% Sig_chol
  
  for(tt in 2:nirf){
    cc_dy <- B_comp %*% cc_dy
    virf[tt, ,] <- t(msel) %*% cc_dy %*% msel %*% Sig_chol
  }

  VAR$Sig_chol <- Sig_chol
  VAR$cholIRF <- virf

  return(VAR)
}

chol_VAR <- function(VAR){
  h <- VAR$h
  
  # varbc <- c()
  # varbc$Ny <- VAR$Ny
  # varbc$p <- VAR$p
  # varbc$Perm <- VAR$Perm
  # varbc$irhor <- VAR$irhor
  # varbc$normalize <- VAR$normalize
  # varbc$bet <- VAR$bet
  # varbc$Sigma <- VAR$Sigma ##
  # varbc$res <- VAR$res
  # varbc
  var_chol_bc <- irfvar(VAR)
  return(var_chol_bc)
}

wild_bootstrap <- function(VAR, nboot, clevel) {
  # Detrend the VAR residuals
  res <- detrend(VAR$res, tt = "constant")
  
  VARbs <- list()
  VARbs$irf <- array(0, dim = c(VAR$irhor, VAR$Ny, VAR$Ny, nboot))
  VARbc <- list()
  
  for (jj in 1:nboot) {
    # Resample full cross-section of residuals within a time period
    Tmax <- max(VAR$Ti)
    cal <- seq(min(VAR$Year), max(VAR$Year))
    rr <- sample(c(-1, 1), size = Tmax, replace = TRUE)
    resb <- matrix(0, nrow = VAR$obs, ncol = VAR$Ny)
    
    for (nn in 1:VAR$obs) {
      year_idx <- which(cal == VAR$Year[nn])
      resb[nn, ] <- rr[year_idx] * res[nn, ]
    }
    
    t <- 1
    Yb <- matrix(0, nrow = VAR$obs, ncol = VAR$Ny)
    Xb <- matrix(0, nrow = VAR$obs, ncol = VAR$Ny * VAR$p)
    
    for (i in 1:max(VAR$ii)) {
      Xb[t, ] <- VAR$X[t, 1:(VAR$Ny * VAR$p)]
      for (j in t:(t + VAR$Ti[i] - 1)) {
        Yb[j, ] <- Xb[j, ] %*% VAR$bet_bc[1:(VAR$Ny * VAR$p), , drop = FALSE] +
          c(VAR$D[j, 1:(VAR$N - 1)], 1) %*% VAR$bet_bc[(VAR$Ny * VAR$p + 1):nrow(VAR$bet_bc), , drop = FALSE] +
          resb[j, ]
        if (j < t + VAR$Ti[i] - 1) {
          Xb[j + 1, ] <- c(Yb[j, ], Xb[j, 1:((VAR$p - 1) * VAR$Ny)])
        }
      }
      t <- t + VAR$Ti[i]
    }
    
    # Run VAR with this bootstrapped data
    VARBS <- list()
    VARBS$I_m <- VAR$I_m
    VARBS$p <- VAR$p
    VARBS$Ny <- VAR$Ny
    VARBS$N <- VAR$N
    VARBS$irhor <- VAR$irhor
    VARBS$Perm <- VAR$Perm
    VARBS$normalize <- VAR$normalize
    VARBS$Y <- Yb
    VARBS$X <- cbind(Xb, VAR$D[, 1:(VAR$N - 1)])
    
    VARBS$bet <- qr.solve(cbind(VARBS$X, rep(1, VAR$obs)), VARBS$Y)
    
    # Bias-correct the estimate
    VARBS$bet[1:(VAR$Ny * VAR$p), ] <- VARBS$bet[1:(VAR$Ny * VAR$p), ] - VAR$bc_diff
    VARBS$res <- VARBS$Y - cbind(VARBS$X, rep(1, VAR$obs)) %*% VARBS$bet
    VARBS$Sigma <- t(VARBS$res) %*% VARBS$res / (VAR$obs - VAR$Ny * VAR$p - VAR$N)
    
    # Compute Cholesky impulse response function
    VARBS <- irfvar(VARBS)
    
    VARbs$irf[, , , jj] <- VARBS$cholIRF
  }
  
  # IRF confidence intervals for bias-corrected VAR
  VARbc$CholirsH1 <- apply(VARbs$irf, c(1, 2, 3), function(x) quantile(x, (1 - clevel[1] / 100) / 2))
  VARbc$CholirsL1 <- apply(VARbs$irf, c(1, 2, 3), function(x) quantile(x, 1 - (1 - clevel[1] / 100) / 2))
  VARbc$CholirsH2 <- apply(VARbs$irf, c(1, 2, 3), function(x) quantile(x, (1 - clevel[2] / 100) / 2))
  VARbc$CholirsL2 <- apply(VARbs$irf, c(1, 2, 3), function(x) quantile(x, 1 - (1 - clevel[2] / 100) / 2))
  VARbc$Cholirsmean <- apply(VARbs$irf, c(1, 2), mean)
  VARbc$Xb <- Xb
  return(VARbc)
}


plot_irf <- function(vec_irf, vec_hci, vec_lci){
  len <- length(vec_irf)
  g <- data.frame(t = 0:(len-1), irf = vec_irf, hci = vec_hci, lci = vec_lci) %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(aes(x = t, y = irf), shape = 8, col = "red") +
    ggplot2::geom_line(aes(x = t, y = irf), col = "red") +
    ggplot2::geom_line(aes(x = t, y = hci), col = "red", linetype = "dashed") +
    ggplot2::geom_line(aes(x = t, y = lci), col = "red", linetype = "dashed")
}

regress_lagged_y <- function(data, i){
  data <- data %>% 
    group_by(c) %>% 
    mutate(lnGDP_d3_lagged = dplyr::lead(lnGDP_d3, i, order_by = c))
  result <- fixest::feols(lnGDP_d3_lagged ~ D3HHD_GDP + D3NFD_GDP | c,
                          data = data)
  return(result)
}
