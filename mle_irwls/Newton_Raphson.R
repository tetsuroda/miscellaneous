# ##########################################
# Newton-Raphson Iteration for poisson MLE
# ##########################################

# Function for Newton-Raphson
poisson_MLE <-
  
  function(mu_int,y, Iter) {
    
    mu <- matrix(NA, nrow = Iter+1, ncol=1)
    mu[1] <- mu_int
    dl <- matrix(NA, nrow = Iter+1, ncol=1)
    ddl <- matrix(NA, nrow = Iter+1, ncol=1)
    
    for (i in 1:Iter){
      
      dl[i] <- sum(y)/mu[i] - length(y)
      ddl[i] <- (-mu[i]^-2)*sum(y)
      
      mu[i+1] <- mu[i] - dl[i] / ddl[i]
      
    }
    return(data.frame(mu=mu,score=dl,sec_dl=ddl))
  }


N <- 100
y <- rpois(n=N,lambda=5)

# mu and other results
iter_res <- poisson_MLE(1,y,20)

# Visualise
library(ggplot2)
ggplot(iter_res, aes(x=1:nrow(iter_res),y=mu)) + geom_line()
ggplot(iter_res, aes(x=1:nrow(iter_res),y=score)) + geom_line()
ggplot(iter_res, aes(x=1:nrow(iter_res),y=sec_dl)) + geom_line()


# #####################################################
# Newton-Raphson Iteration for poisson regression models
# #####################################################

# Toy data
# -------------------------------------------
N <- 100

x1 <- runif(n=N, min = -5, max = 5)
x2 <- rnorm(n=N, mean= 1, sd = 1)

X <- matrix(c(rep(1,N),x1, x2), nrow = N)
b_true <- c(1,1,1)
xb <- X%*%b_true

library(tidyverse)
# log(u) = 1 + x1 + x2 -> u = exp(1+x1+x2)
y <- data.frame(x=exp(xb)) %>%
  pmap_dbl(~rpois(n=1, lambda=.x))

df <- data.frame(x1,x2,y)

# Visualisation
# ---------------------------------------------
df %>% ggplot(aes(x=x1,y=y)) + geom_point()
df %>% ggplot(aes(x=x2,y=y)) + geom_point()
df %>% ggplot(aes(x=y)) + geom_histogram()


# Function
# ---------------------------------------------
irwls_fun <-
  function(X,b_int,y,Iter) {
    
    N <- nrow(X)
    l <- matrix(NA, ncol=1, nrow=Iter+1)
    b <- matrix(NA, ncol=3, nrow=Iter+1)
    b[1,] <- b_int
    # log likelihood
    l[1] <- t(y)%*%X%*%b[1,] - sum(exp(X%*%b[1,]))
    
    for (i in 1:Iter){
      # Weights
      # Model: log(u) = b0 + b1*x1 + b2*x2
      # u = exp(b[i,1]+b[i,2]*x1+b[i,3]*x2)
      # var(y) = u
      # d u_i / d eta_i = exp(eta_i) = u_i
      u <- exp(X%*%b[i,]) %>% as.vector()
      W <-diag((1/u)*u^2, ncol = N)

      # Information matrix
      t(X)%*%W%*%X
      
      # z
      z <- X%*%b[i,] + (y-u)*u^-1

      # b
      b[i+1,] <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z %>% t()
      
      # log likelihood
      l[i+1] <- t(y)%*%X%*%b[i+1,] - sum(exp(X%*%b[i+1,]))
    }
    
    res <- cbind(seq(0,Iter,1),b,l)
    colnames(res) <- c("Iter","b0","b1","b2","loglik")
    return(res %>% data.frame())
    
  }


# Run the function
# -----------------------------------------------------

# The length until the convergence varies by initial values
res <- irwls_fun(X = X, b_int = c(3,3,3), y = df$y, Iter = 30)

res %>% ggplot(aes(x=Iter, y=b0)) + geom_line()
res %>% ggplot(aes(x=Iter, y=b1)) + geom_line()
res%>% ggplot(aes(x=Iter, y=b2)) + geom_line()

# With the initial values, they also converge around 1, the true parameter value.
res <- irwls_fun(X = X, b_int = c(-10,-10,-10), y = df$y,Iter = 100)
res %>% ggplot(aes(x=Iter, y=b0)) + geom_line()
res %>% ggplot(aes(x=Iter, y=b1)) + geom_line()
res %>% ggplot(aes(x=Iter, y=b2)) + geom_line()