calculate_theta_coefs <- function(w, X, Y, C, D) {
  return(Inverse(w * t(X) %*% X + (1 - w) * t(Y) %*% Y, tol = 1e-14) %*% 
         (w * t(X) %*% C + (1 - w) * t(Y) %*% D))
}

calculate_delta <- function(S0, sigma, Tt, w, K = 1, pol_degree = 7){
  
  rand_norm_T <- rnorm(length(S0))
  S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
  
  zeroes <- c(rep(0, length(S0)))
  
  ones <- c(rep(1, length(S0)))
  
  Y <- zeroes %>% 
    cbind(ones) %>% 
    cbind(poly(S0, pol_degree-1, raw=TRUE)) %>% 
    as.matrix()
  
  for(i in 1:(pol_degree+1)){
    Y[,i] <- Y[,i] * (i-1)
  }
  
  C <- S_T
  
  X <- ones %>%
    cbind(poly(S0, pol_degree, raw=TRUE)) %>%
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

est_derr_pricing_function <- function(coefs, S0, pol_degree = 7){
  S0_6 <- poly(S0, (pol_degree-1), raw = TRUE)
  return((S0_6 %*% (c(2:pol_degree) * coefs[2:pol_degree])) + coefs[1])
}


calculate_delta_price <- function(S0, sigma, Tt, K = NULL, w = NULL, pol_degree = 7){
  rand_norm_T <- rnorm(length(S0))
  S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
  
  call_prices <- sapply(S_T, FUN = function(x)max(x - K, 0))
  
  # LM
  call_simulations <- cbind(S0, call_prices) %>% 
    as_data_frame()
  
  est_pricing_func <- lm(call_prices ~ poly(S0, pol_degree, raw=TRUE), data = call_simulations)
  
  coefs <- est_pricing_func$coefficients %>% as_tibble() %>% mutate(value = ifelse(is.na(value), 0, value)) %>% as.matrix()

  derr_pricing_func_val <- est_derr_pricing_function(coefs = coefs[-1],
                                                     S0 = S0, 
                                                     pol_degree = (pol_degree - 1))
  return(list(delta = derr_pricing_func_val))
}

true_delta <- function(S0, K, sigma, Tt, w, pol_degree){
  return(list(delta = pnorm((S0-K)/(sigma * sqrt(Tt)))))
}

calculate_hedge_error <- function(dt, Tt, num_rep, K, sigma, St, delta_func, S0, w = NULL, pol_degree = 7){
  S0 <- S0
  for(t in seq(dt, Tt, dt)){
    St <- St + sigma * sqrt(dt) * rnorm(num_rep)
    Vt <- at * St + bt
    at <- delta_func(S0 = S0, Tt = (1-t), K = K, sigma = sigma, w = w, pol_degree = 7)$delta
    bt <- Vt - at * St
  }
  
  err <- sd((Vt - sapply(St, FUN = function(x)max(x - K, 0))))
  return(err)
}
