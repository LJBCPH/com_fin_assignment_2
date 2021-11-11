calculate_theta_coefs <- function(w, X, Y, C, D) {
  return(Inverse(w * t(X) %*% X + (1 - w) * t(Y) %*% Y, tol = 1e-14) %*% 
         (w * t(X) %*% C + (1 - w) * t(Y) %*% D))
}

calculate_delta <- function(S0, sigma, Tt, w, K = 1, pol_degree = 7, St = NULL, coefs = NULL, BS = FALSE, LRM = FALSE, digital = FALSE, MM = FALSE, eps = 10, setseed = F){
  if(is.null(coefs)){
    if(setseed == T){
      set.seed(1)
    }
    rand_norm_T <- rnorm(length(S0))
    if(BS == FALSE){
      S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
      
      D <- S_T %>% 
        as_tibble() %>% 
        transmute(D_val = ifelse(value >= K, 1, 0)) %>% 
        pull()
    } else {
      S_T <- S0 * exp(-1/2 * sigma^2*Tt + sigma * sqrt(Tt) * rand_norm_T)
      if(digital == TRUE){
        payoff <- ifelse(S_T >= K, 1, 0)
      } else {
        payoff <- (S_T - K) * ifelse(S_T >= K, 1, 0)
      }
      if(LRM == FALSE) {
        D <- cbind(S_T, S0) %>% 
          as_tibble() %>% 
          transmute(D_val = S_T * ifelse(S_T >= K, 1, 0) / S0) %>% 
          pull()
      } else {
        if(MM == FALSE){
          D <- cbind(S_T, S0) %>% 
            as_tibble() %>% 
            transmute(D_val = payoff * rand_norm_T / (sigma * sqrt(Tt) * S0)) %>% 
            pull()
        } else {
          D <- cbind(S_T, S0) %>% 
            as_tibble() %>% 
            transmute(D_val = (1/(2 * eps) * (abs(S_T-K)<eps) * S_T/S0 + 
                                 ((S_T>K)-(pmin(1, pmax(0, S_T - K + eps)/(2*eps)))) *
                                 rand_norm_T/(S_0*sqrt(Tt)*sigma))) %>% 
            pull()
        }
      }
    }
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
  
    coefs <- calculate_theta_coefs(w, X, Y, C, D)
  }
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
                D = D,
                coefs = coefs))
  }
}

est_derr_pricing_function <- function(coefs, S0, pol_degree = 7){
  S0_6 <- poly(S0, (pol_degree-1), raw = TRUE)
  return((S0_6 %*% (c(2:pol_degree) * coefs[2:pol_degree])) + coefs[1])
}

calculate_delta_price <- function(S0, sigma, Tt, K = NULL, w = NULL, pol_degree = 7, St = NULL, 
                                  coefs = NULL, BS = FALSE, LRM = FALSE, digital = F, MM = FALSE, eps = NULL){
  #set.seed(1)
  rand_norm_T <- rnorm(length(S0))
  #S0 <- 1 + sigma * sqrt(Tt) * rand_norm_T
  if(BS == FALSE){
    S_T <- S0 + sigma * sqrt(Tt) * rand_norm_T
  } else {
    S_T <- S0 * exp((1/2 * sigma^2)*Tt + sigma * sqrt(Tt) * rand_norm_T)
  }
  if(digital == T){
    call_prices <- sapply(S_T, FUN = function(x)ifelse(x >= K, 1, 0))
  } else {
    call_prices <- sapply(S_T, FUN = function(x)max(x - K, 0))
  }
  
  # LM
  call_simulations <- cbind(S0, call_prices) %>% 
    as_data_frame()
  
  est_pricing_func <- lm(call_prices ~ poly(S0, pol_degree, raw=TRUE), data = call_simulations)
  
  coefs <- est_pricing_func$coefficients %>% as_tibble() %>% mutate(value = ifelse(is.na(value), 0, value)) %>% as.matrix()

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

true_delta <- function(S0 = NULL, K, sigma, Tt, w = NULL, pol_degree = NULL, St, 
                       coefs = NULL, BS = FALSE, LRM = FALSE, digital = FALSE, eps = NULL, MM = FALSE){
  if(digital == FALSE){
    if(BS == FALSE){
      return(list(delta = pnorm((St-K)/(sigma * sqrt(Tt)))))
    } else {
      return(list(delta = pnorm((log(St / K) + (1/2 * sigma^2) * Tt) / (sigma * sqrt(Tt)))))
    }
  } else {
    d1 <- (log(St / K) + 1/2 * sigma^2 * Tt) / (sigma * sqrt(Tt))
    return(list(delta = (dnorm(d1 - sigma * sqrt(Tt))) / (St * sqrt(Tt) * sigma)))
  }
}

ZeroRateBachelierCall <- function(St, Tt, K, sigma){
  d<-(St-K)/(sigma*sqrt(Tt))
  CallDelta<-pnorm(d,0,1)
  CallPrice<-(St-K)*CallDelta+sigma*sqrt(Tt)*dnorm(d,0,1)
  return(list(Price=CallPrice, delta=CallDelta))
}

digital_call <- function(St, Tt, K, sigma){
  d1 = (log(St / K) + 1/2 * sigma^2 * Tt) / (sigma * sqrt(Tt))
  d2 = d1 - sigma * sqrt(Tt)
  return(list(Price = pnorm(d2)))
}

calculate_hedge_error <- function(dt = 1/52, Tt, num_rep, K, sigma, St, S0, delta_func, w = NULL, pol_degree = 7, seed = 3, 
                                  call_func = ZeroRateBachelierCall, coefs = NULL, BS = FALSE, LRM = FALSE, digital = FALSE, MM = FALSE,
                                  eps = 1){
  set.seed(seed)
  # Initializing values
  initial_price <- call_func(St = St, Tt = Tt, K = K, sigma = sigma)$Price
  Vt <- initial_price
  at <- true_delta(S0 = S0, Tt = Tt, K = K, sigma = sigma, w = w, pol_degree = pol_degree, 
                   St = St, BS = BS, LRM = LRM, digital = digital, MM = MM, eps = eps)$delta
  bt <- Vt-at*St
  rand_norm_0 <- rnorm(num_sim)
  
  if(BS == FALSE){
    S0 <- K + d * sigma * sqrt(Tt) * rand_norm_0
  }else{
    S_0 <- K * exp(-1/2 * sigma^2 * Tt + sigma * sqrt(Tt) * rand_norm_0)
  }
  
  ### LOOPING THROUGH ALL TIME TO MATURITIES AND UPDATING VALUE OF HEDGEPORTFOLIO
  for(t in seq(dt, Tt-dt, dt)){
    if(BS == FALSE){
      St <- St + sigma * sqrt(dt) * rnorm(num_rep)
    } else {
      St <- St * exp(-1/2 * sigma^2 * dt + sigma * sqrt(dt) * rnorm(num_rep))
    }
    Vt <- at * St + bt
    at <- delta_func(S0 = S0, Tt = 1-t, K = K, sigma = sigma, w = w, pol_degree = pol_degree, St = St, 
                     coefs = coefs, BS = BS, LRM = LRM, digital = digital, MM = MM, eps = eps)$delta
    bt <- Vt - at * St
  }
  
  if(BS == FALSE) {
    St <- St + sigma * sqrt(dt) * rnorm(num_rep)
  } else {
    St <- St * exp(-1/2 * sigma^2 * dt + sigma * sqrt(dt) * rnorm(num_rep))
  }
  Vt <- at * St + bt
  if(digital == TRUE){
    err <- sd(Vt - sapply(St, FUN = function(x)ifelse(x >= K, 1, 0)))
  } else {
    err <- sd(Vt - sapply(St, FUN = function(x)max(x - K, 0)))
  }
  return(err)
}

BS_func <- function(St, K, sigma, Tt){
  d1 <- (log(St / K) + (sigma^2/2)*Tt) / (sigma * sqrt(Tt))
  d2 <- d1 - sigma * sqrt(Tt)
  value <- St * pnorm(d1) - K * pnorm(d2)
  return(list(Price = value))
}
