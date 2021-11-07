calculate_theta_coefs <- function(w, X, Y, C, D) {
  return(Inverse(w * t(X) %*% X + (1 - w) * t(Y) %*% Y, tol = 1e-14) %*% 
         (w * t(X) %*% C + (1 - w) * t(Y) %*% D))
}

calculate_delta <- function(S0, sigma, Tt, w, K = 1, pol_degree = 7, St = NULL){
  
  rand_norm_T <- rnorm(length(S0))
  #S0 <- 1 + sigma * sqrt(Tt) * rand_norm_T
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

  coefs <- calculate_theta_coefs(w, X, Y, C, D)
  if(missing(St)){
    return(list(delta = Y %*% coefs, 
                D = D))
  }else{
    zeroes <- c(rep(0, length(St)))
    
    ones <- c(rep(1, length(St)))
    
    Y <- zeroes %>% 
      cbind(ones) %>% 
      cbind(poly(St, pol_degree-1, raw=TRUE)) %>% 
      as.matrix()
    
    for(i in 1:(pol_degree+1)){
      Y[,i] <- Y[,i] * (i-1)
    }
    
    return(list(delta = Y %*% coefs, 
                D = D))
  }
}

est_derr_pricing_function <- function(coefs, S0, pol_degree = 7){
  S0_6 <- poly(S0, (pol_degree-1), raw = TRUE)
  return((S0_6 %*% (c(2:pol_degree) * coefs[2:pol_degree])) + coefs[1])
}

calculate_delta_price <- function(S0, sigma, Tt, K = NULL, w = NULL, pol_degree = 7, St = NULL){
  rand_norm_T <- rnorm(length(S0))
  #S0 <- 1 + sigma * sqrt(Tt) * rand_norm_T
  S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
  
  call_prices <- sapply(S_T, FUN = function(x)max(x - K, 0))
  
  # LM
  call_simulations <- cbind(S0, call_prices) %>% 
    as_data_frame()
  
  est_pricing_func <- lm(call_prices ~ poly(S0, pol_degree, raw=TRUE), data = call_simulations)
  
  coefs <- est_pricing_func$coefficients %>% as_tibble() %>% mutate(value = ifelse(is.na(value), 0, value)) %>% as.matrix()
  # if(Tt < 0.15){
  # cat("polynm: ", pol_degree, " t: ",  Tt, " coefs: ", coefs, "\n")
  # }
  if(missing(St)){
    derr_pricing_func_val <- est_derr_pricing_function(coefs = coefs[-1],
                                                       S0 = S0, 
                                                       pol_degree = (pol_degree))
  } else {
    derr_pricing_func_val <- est_derr_pricing_function(coefs = coefs[-1],
                                                       S0 = St, 
                                                       pol_degree = (pol_degree))
  }
  return(list(delta = derr_pricing_func_val))
}

true_delta <- function(S0, K, sigma, Tt, w, pol_degree, St){
  return(list(delta = pnorm((St-K)/(sigma * sqrt(Tt)))))
}

ZeroRateBachelierCall<-function(S,T,K,vol){
  d<-(S-K)/(vol*sqrt(T))
  CallDelta<-pnorm(d,0,1)
  CallPrice<-(S-K)*CallDelta+vol*sqrt(T)*dnorm(d,0,1)
  return(list(Price=CallPrice, Delta=CallDelta))
}

calculate_hedge_error <- function(dt, Tt, num_rep, K, sigma, St, S0, delta_func, w = NULL, pol_degree = 7, seed = 3){
  set.seed(seed)
  initial_price <- ZeroRateBachelierCall(St, Tt, K, sigma)$Price
  Vt <- initial_price
  at <- true_delta(S0 = S0, Tt = Tt, K = K, sigma = sigma, w = w, pol_degree = pol_degree, St = St)$delta
  bt <- Vt-at*St
  
  for(t in seq(dt, Tt-dt, dt)){
    rand_norm_0 <- rnorm(num_sim)
    S_0 <- K + d * sigma * sqrt(Tt) * rand_norm_0
    St <- St + sigma * sqrt(dt) * rnorm(num_rep)
    Vt <- at * St + bt
    at <- delta_func(S0 = S0, Tt = Tt-t, K = K, sigma = sigma, w = w, pol_degree = pol_degree, St = St)$delta
    bt <- Vt - at * St
  }
  St <- St + sigma * sqrt(dt) * rnorm(num_rep)
  Vt <- at * St + bt
  err <- sd(Vt - sapply(St, FUN = function(x)max(x - K, 0)))
  return(err)
}

for (i in 2:(Nhedge-1)) {
  S<-S+vol*sqrt(dt)*rnorm(Npaths)
  tau<-T-dt*(i-1)
  dummy<-ZeroRateBachelierCall(S,tau,K,vol)
  if (Poly){
    for (k in 1:Npaths)dummy$Delta[k]<-CoefMatrix[2:M,i]%*%(Powers[2:M]*S[k]^(Powers[2:M]-1))
    
  }
  Vpf<-a*S+b
  a<-dummy$Delta
  b<-Vpf-a*S
}
