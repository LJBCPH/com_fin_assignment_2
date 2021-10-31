calculate_theta_coefs <- function(w, X, Y, C, D) {
  return(Inverse(w * t(X) %*% X + (1 - w) * t(Y) %*% Y, tol = 1e-14) %*% 
         (w * t(X) %*% C + (1 - w) * t(Y) %*% D))
}

calculate_delta <- function(S0, sigma, Tt, w){
  
  rand_norm_T <- rnorm(length(S0))
  S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
  
  zeroes <- c(rep(0, length(S0)))
  
  ones <- c(rep(1, length(S0)))
  
  Y <- zeroes %>% 
    cbind(ones) %>% 
    cbind(poly(S0, 6, raw=TRUE)) %>% 
    as.matrix()
  
  for(i in 1:8){
    Y[,i] <- Y[,i] * (i-1)
  }
  
  C <- S_T
  
  X <- ones %>%
    cbind(poly(S0, 7, raw=TRUE)) %>%
    as.matrix()
  
  D <- S_T %>% 
    as_tibble() %>% 
    transmute(D_val = ifelse(value >= K, 1, 0)) %>% 
    pull()
  # cat("dim X: ", dim(X), 
  #     "\n dim Y: ", dim(Y),
  #     "\n dim C: ", length(C))
  return(list(delta = Y %*% calculate_theta_coefs(w, X, Y, C, D), 
              D = D))
}

est_derr_pricing_function <- function(coefs, S0){
  S0_6 <- poly(S0, 6, raw = TRUE)
  return((S0_6 %*% (c(2:7) * coefs[2:7])) + coefs[1])
}

calculate_delta_price <- function(S0, sigma, Tt){
  rand_norm_T <- rnorm(length(S0))
  S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
  
  call_prices <- sapply(S_T, FUN = function(x)max(x - K, 0))
  
  # LM
  call_simulations <- cbind(S0, call_prices) %>% 
    as_data_frame()
  
  est_pricing_func <- lm(call_prices ~ poly(S0, 7, raw=TRUE), data = call_simulations)
  
  coefs <- est_pricing_func$coefficients %>% as_tibble() %>% mutate(value = ifelse(is.na(value), 0, value)) %>% as.matrix()
  
  derr_pricing_func_val <- est_derr_pricing_function(coefs[-1],
                                                     S0)
  return(derr_pricing_func_val)
}

true_delta <- function(St, K, sigma, t){
  return(pnorm((St-K)/(sigma * sqrt(t))))
}

calculate_hedge_error <- function(dt, Tt, num_rep, K, sigma, St, delta_func, ...){
  for(t in seq(dt, Tt, dt)){
    St <- St + sigma * sqrt(dt) * rnorm(num_rep)
    Vt <- at * St + bt
    at <- delta_func(St = St, t = t, ...)
    bt <- Vt - at * St
  }
  
  err <- sd((Vt - sapply(St, FUN = function(x)max(x - K, 0))))
  return(err)
}
