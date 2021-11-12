################## ROLFS CODE ###################
tridagsoln<-function(a,b,c,r) {	
  
  n<-length(b)
  gam<-u<-tridagsoln<-1:n
  bet<-b[1]
  u[1]<-r[1]/bet
  
  for (j in 2:n) {
    gam[j]<-c[j-1]/bet
    bet<-b[j]-a[j]*gam[j]
    u[j]<-(r[j]-a[j]*u[j-1])/bet
  }
  
  for(j in (n-1):1){ 
    u[j]<-u[j]- gam[j+1]*u[j+1]
  }
  return(u)
}
################################################

################## Andreasens CODE ###################
tridagsoln<-function(a,b,c,r) {	
  
  dim <- length(b)
  ws <- c(rep(0, dim))
  bet <- 1 / b[1]
  x[1] <- r[1] * bet
  
  for (j in 2:dim){
    ws[j] <- c[j-1] * bet
    bet <- 1 / (b[j] - a[j] * ws[j])
    x[j] <- (x[j] - a[j] * x[j-1]) * bet
  }
  
  for(j in (n-1):2){ 
    x[j] <- x[j] - ws[j+1] * x[j+1]
  }
  return(x)
}
################################################

################# FROM CLAUS MUNK ###############
tri_diag_solver <- function(a, b, c, d, payoff = NULL) {	
  dim <- length(b)
  f <- rep(0, dim)
  for(j in 1:(dim - 1)){
    #Gaussian elimination
    a[j+1] <- a[j+1] / b[j]
    b[j+1] <- b[j+1] - a[j+1] * c[j]
    d[j+1] <- d[j+1] - a[j+1] * d[j] #Forward substitution
  }
  if(!is.null(payoff)){
    f[dim] <- d[dim] / b[dim]
  } else {
    f[dim] <- max(d[dim] / b[dim], payoff[dim])
  }
  for(i in (dim-1):1){
    f[i] <- max((d[i] - c[i] * f[i+1]) / b[i], payoff[i]) #Backward substitution
  }
  return(f)
}
################################################

###################### CRANK NICOLSON EU####################
crank_nicolson_eu <- function(r, Tt, x_min, x_max, time_steps, x_steps, sigma, K = 40, spots = c(36, 38, 40, 42, 44), gamma = 1, dupire = F, selected_times = c(1, 0)){
  dt <- Tt / time_steps
  dx <- (x_max - x_min) / x_steps
  grid_time <- seq(0, dt * time_steps, dt)
  grid_x <- seq(0, x_min + dx * x_steps, dx)
  rhs_step <- c(rep(0, x_steps - 1))
  
  mu <- r * grid_x
  sigma_2 <- signs * sigma^2 * grid_x^(2*gamma)
  if(dupire == FALSE) {
    signs <- 1
  } else {
    signs <- -1
  }
  
    a_lhs <- 1/4 * dt * (sigma_2 / (dx^2) - mu / dx)
    b_lhs <- -1 - 1/2 * dt * (sigma_2 / (dx^2) + r)
    c_lhs <- 1/4 * dt * (sigma_2 / (dx^2) + mu / dx)
    b_star <- -1 + 1/2 * dt * (sigma_2 / (dx^2) + r)
  if(dupire == TRUE){
    temp <- b_lhs
    b_lhs <- b_star
    b_star <- temp
    temp <- K
    K <- spots
    spots <- temp
    rm(temp)
  }
    
    # Defining matrix form like eq. 15, using the right hand side from eq. 19
    rhs_m <- diag(b_star[-c(1, (x_steps + 1))]) # diagonal
    rhs_m[row(rhs_m) - col(rhs_m) == 1] <- signs * -a_lhs[3:(x_steps)] # lower diagona l
    rhs_m[col(rhs_m) - row(rhs_m) == 1] <- signs * -c_lhs[2:(x_steps - 1)] # Upper diagonal
    
    # Defining our f-vector from eq. 19 in Munks Note
    f <- matrix(0, nrow = x_steps + 1, ncol = time_steps + 1)
    if(dupire == FALSE){
      f[1,] <- exp(-r*(Tt - grid_time)) * (K - x_min)
      f[, time_steps + 1] <- sapply(grid_x, FUN = function(x)max(0, K - x))
    } else {
      f[1,] <- exp(r * (grid_time)) * K - x_min
      f[,1] <- sapply(grid_x, FUN = function(x)max(0, K - x))
    }
    
    # Removing first and last input to create the shown matrix on the left hand side (eq. 19)
    a_lhs_mod <- signs * a_lhs[-c(1, x_steps + 1)]
    b_lhs_mod <- signs * b_lhs[-c(1, x_steps + 1)]
    c_lhs_mod <- signs * c_lhs[-c(1, x_steps + 1)]

  for(i in time_steps:1){
    rhs_step[1] <- signs * b_lhs[1] * f[1, i] + signs * c_lhs[1] * f[2, i]
    rhs_step[x_steps - 1] <- signs * a_lhs[x_steps + 1] * f[x_steps, i] + signs * b_lhs[x_steps + 1] * f[x_steps + 1, i]
    rhs <- rhs_m %*% f[-c(1, x_steps + 1), i + signs * 1] + rhs_step
    
    f[-c(1, x_steps + 1), i] <- tri_diag_solver(a_lhs_mod, b_lhs_mod, c_lhs_mod, rhs)
  }
  if(dupire == T){
    for (i in 2:(time_steps + 1)){
      f[,i] <- exp(-r * grid_time[i]) * f[,i]
    }
  }  
  approx_spots <- grid_time %>% rbind(grid_x %>% as_tibble() %>% bind_cols(f %>% as_tibble()))
  approx_spots <- sapply(spots, FUN = function(x) approx_spots[which(abs(x - approx_spots$value) == min(abs(x - approx_spots$value))),1])
  approx_spots <- approx_spots %>% unlist() %>% as.data.frame() %>% pull()
  
  approx_times <- sapply(selected_times, FUN = function(x) grid_time[which(abs(x - grid_time) == min(abs(x - grid_time)))])
  rownames(f) <- grid_x; colnames(f) <- grid_time
  res <- f[which(rownames(f) %in% approx_spots), colnames(f) %in% approx_times]

  return(res)
}

########### CRANK NICOLSON US ######################
###################### CRANK NICOLSON EU####################
crank_nicolson_us <- function(r, Tt, x_min, x_max, time_steps, x_steps, sigma, K = 40, spots = c(36, 38, 40, 42, 44), selected_times = c(1, 0)){
  dt <- Tt / time_steps
  dx <- (x_max - x_min) / x_steps
  grid_time <- seq(0, dt * time_steps, dt)
  grid_x <- seq(0, x_min + dx * x_steps, dx)
  rhs_step <- c(rep(0, x_steps - 1))
  
  mu <- r * grid_x
  sigma_2 <- sigma^2 * grid_x^2
  
  a_lhs <- 1/4 * dt * (sigma_2 / (dx^2) - mu / dx)
  b_lhs <- -1 - 1/2 * dt * (sigma_2 / (dx^2) + r)
  c_lhs <- 1/4 * dt * (sigma_2 / (dx^2) + mu / dx)
  b_star <- -1 + 1/2 * dt * (sigma_2 / (dx^2) + r)
  
  # Defining matrix form like eq. 15, using the right hand side from eq. 19
  rhs_m <- diag(b_star[-c(1, (x_steps + 1))]) # diagonal
  rhs_m[row(rhs_m) - col(rhs_m) == 1] <- -a_lhs[3:(x_steps)] # lower diagona l
  rhs_m[col(rhs_m) - row(rhs_m) == 1] <- -c_lhs[2:(x_steps - 1)] # Upper diagonal
  
  # Defining our f-vector from eq. 19 in Munks Note
  f <- matrix(0, nrow = x_steps + 1, ncol = time_steps + 1)
  f[1,] <- exp(-r*(Tt - grid_time)) * (K - x_min)
  f[, time_steps + 1] <- sapply(grid_x, FUN = function(x)max(0, K - x))
  
  # Removing first and last input to create the shown matrix on the left hand side (eq. 19)
  a_lhs_mod <- a_lhs[-c(1, x_steps + 1)]
  b_lhs_mod <- b_lhs[-c(1, x_steps + 1)]
  c_lhs_mod <- c_lhs[-c(1, x_steps + 1)]
  
  payoff <- pmax(K - grid_x[-c(1, (x_steps + 1))], 0)
  for(i in time_steps:1){
    f[1, i] <- max(exp(-r*(Tt - i * dt)) * f[1, i + 1], K - x_min)
    rhs_step[1] <- b_lhs[1] * f[1, i] + c_lhs[1] * f[2, i]
    rhs_step[x_steps - 1] <- a_lhs[x_steps + 1] * f[x_steps, i] + b_lhs[x_steps + 1] * f[x_steps + 1, i]
    rhs <- rhs_m %*% f[-c(1, x_steps + 1), i + 1] + rhs_step
    
    f[-c(1, x_steps + 1), i] <- tri_diag_solver(a_lhs_mod, b_lhs_mod, c_lhs_mod, rhs, payoff)
  }
  rownames(f) <- grid_x; colnames(f) <- grid_time
  
  approx_spots <- grid_time %>% rbind(grid_x %>% as_tibble() %>% bind_cols(f %>% as_tibble()))
  approx_spots <- sapply(spots, FUN = function(x) approx_spots[which(abs(x - approx_spots$value) == min(abs(x - approx_spots$value))),1])
  approx_spots <- approx_spots %>% unlist() %>% as.data.frame() %>% pull()
  
  approx_times <- sapply(selected_times, FUN = function(x) grid_time[which(abs(x - grid_time) == min(abs(x - grid_time)))])
  
  res <- f[which(rownames(f) %in% approx_spots), colnames(f) %in% approx_times]
  return(res)
  
}
##################################################

