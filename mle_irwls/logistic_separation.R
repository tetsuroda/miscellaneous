library(tidyverse)

# Prepare a dataset
N <- 10
x1 <- rnorm(n=N, mean = 1, sd = 3)
X <- matrix(c(rep(1,N),x1), nrow = N, ncol = 2)
b_true <- c(1,2)
xb <- X%*%b_true

# log(u) = 1 + x1 + x2 -> u = exp(1+x1+x2)
y <- data.frame(x=(exp(xb)/(1+exp(xb)))) %>%
  pmap_dbl(~rbinom(n = 1, size = 1, prob =.x))

df <- data.frame(x1,y)

# Fit a logistic model using glm()
fit <- glm(y ~ x1, data = df, family = "binomial")
coef(fit)


# Define my own logistic function
irwls_logistic <-
  function(X,b_int,y,Iter) {
    
    N <- nrow(X)
    rep_num <- matrix(NA, ncol=1, nrow=Iter+1)
    l <- matrix(NA, ncol=1, nrow=Iter+1)
    b <- matrix(NA, ncol=2, nrow=Iter+1)
    p_x1_0 <- matrix(NA, ncol=1, nrow=Iter+1)
    p_x1_1 <- matrix(NA, ncol=1, nrow=Iter+1)
    
    rep_num[1] <- 1
    
    b[1,] <- b_int
     
    # probability of y being 1
    p <- exp(X%*%b[1,])/(1+exp(X%*%b[1,]))
    
    # log likelihood
    l[1] <- t(y)%*%X%*%b[1,] + sum(log(1 - p))
    
    # probability of y when x1 = 0 or 1
    p_x1_0[1] <- exp(c(1,0)%*%b[1,])/(1+exp(c(1,0)%*%b[1,]))
    p_x1_1[1] <- exp(c(1,1)%*%b[1,])/(1+exp(c(1,1)%*%b[1,]))
    
    for (i in 1:Iter){
      # Weights
      # Model: logit(u) = b0 + b1*x1 
      # logit(pi) = exp(b[i,1]+b[i,2]*x1)
      # u = E(y) = pi 
      # var(y) = N*pi*(1-pi)
      # d u_i / d eta_i = exp(eta_i)/(exp(eta_i) + 1)^2
      # d eta_i / d_u で考えると
      # log(p/(1-p))をpで微分することになるので1/p*(1-p)となりシンプル
      u <- p 
      v <- p*(1-p)
      du_deta <- exp(X%*%b[i,])/(exp(X%*%b[i,]) + 1)^2
      W <- diag(as.vector((1/v)*du_deta^2), ncol = N)
      
      # Information matrix
      J <- t(X)%*%W%*%X
      
      # z
      z <- X%*%b[i,] + (y-u)*du_deta^-1
      
      # b
      b[i+1,] <- solve(J)%*%t(X)%*%W%*%z %>% t()
      
      # probability of y
      p <- exp(X%*%b[i+1,])/(1+exp(X%*%b[i+1,]))
      
      # log likelihood
      l[i+1] <- t(y)%*%X%*%b[i+1,] + sum(log(1 - p))
      
      # probability of y when x1 = 0 or 1
      p_x1_0[i+1] <- exp(c(1,0)%*%b[i+1,])/(1+exp(c(1,0)%*%b[i+1,]))
      p_x1_1[i+1] <- exp(c(1,1)%*%b[i+1,])/(1+exp(c(1,1)%*%b[i+1,]))
      
      # Updated Information matrix
      u <- p 
      v <- p*(1-p)
      du_deta <- exp(X%*%b[i+1,])/(exp(X%*%b[i+1,]) + 1)^2
      W <- diag(as.vector((1/v)*du_deta^2), ncol = N)
      J <- t(X)%*%W%*%X
      J_inv <- solve(J)
      
      rep_num[i+1] <- i+1
      
      if(is.infinite(l[i+1]) == T){
        break
      }
      
      if(abs(l[i+1]-l[i]) < 1e-10) {
        break
      }
        
    }
    
    res <- cbind(rep_num,b,l,p_x1_0,p_x1_1)
    colnames(res) <- c("Iter","b0","b1","loglik","prob_of_y_x_0","prob_of_y_x_1")
    
    iter_res_df <- res %>% data.frame() %>% filter(!is.na(Iter))
    
    estimates <- 
      data.frame(
        beta = c(b0 = b[nrow(iter_res_df),1],
                 b1 = b[nrow(iter_res_df),2]
                 ),
        se = c(b0 = sqrt(J_inv[1,1]),
               b1 = sqrt(J_inv[2,2])
               )
        )
    
    return(list(iter_res_df, estimates))
    
  }

# Fit a logistic regression using my own function
res <- irwls_logistic(X = X, b_int = c(0,0), y = df$y,Iter = 100) 
print(res)

# Plots
plot_loglik <- 
  res[[1]] %>% ggplot(aes(x=Iter, y=loglik)) + 
  geom_line()

plot_b0 <- 
  res[[1]] %>% ggplot(aes(x=Iter, y=b0)) +
  geom_line()

plot_b1 <- 
  res[[1]] %>% ggplot(aes(x=Iter, y=b1)) +
  geom_line()

plot_p_x_0 <-
  res[[1]] %>% ggplot(aes(x=Iter, y=prob_of_y_x_0)) + 
  geom_line() +
  coord_cartesian(ylim = c(0,1)) + ylab("Pr[Y=1|x1=0]")

plot_p_x_1 <-
  res[[1]] %>% ggplot(aes(x=Iter, y=prob_of_y_x_1)) + 
  geom_line() +
  coord_cartesian(ylim = c(0,1)) + ylab("Pr[Y=1|x1=1]")


cowplot::plot_grid(plotlist = list(plot_loglik,
                                   plot_b0,
                                   plot_b1,
                                   plot_p_x_0,
                                   plot_p_x_1))


# ################################
# In case of separation
# ################################

N <- 100
x1 <- rep(c(0,1), each = N/2) 
X <- matrix(c(rep(1,N),x1), nrow = N, ncol = 2)
df <- data.frame(x1) %>%
  mutate(y = ifelse(x1==0,0,1))

# Use glm()
fit <- glm(y ~ x1, data = df, family = "binomial")
summary(fit)
data.frame(df,predict(fit, type = "response")) 

# Use my own function
res <- irwls_logistic(X = X, b_int = c(0,0), y = df$y, Iter = 100) 

# The estimates do not perfectly coincide,
# but it seems this is because they differ on when to stop their iterations
# Mine continues iterating longer than glm() does
# Also, even if looking at the values of b0 and b1 close to 
# those obtained by glm(), I didnt get the exact values.
# I suspect the algorithms used are different.
print(rs)

# Plots
plot_loglik <- 
  res[[1]] %>% ggplot(aes(x=Iter, y=loglik)) + 
  geom_line()

plot_b0 <- 
  res[[1]] %>% ggplot(aes(x=Iter, y=b0)) +
  geom_line()

plot_b1 <- 
  res[[1]] %>% ggplot(aes(x=Iter, y=b1)) +
  geom_line()

plot_p_x_0 <-
  res[[1]] %>% ggplot(aes(x=Iter, y=prob_of_y_x_0)) + 
  geom_line() +
  coord_cartesian(ylim = c(0,1)) + ylab("Pr[Y=1|x1=0]")

plot_p_x_1 <-
  res[[1]] %>% ggplot(aes(x=Iter, y=prob_of_y_x_1)) + 
  geom_line() +
  coord_cartesian(ylim = c(0,1)) + ylab("Pr[Y=1|x1=1]")


cowplot::plot_grid(plotlist = list(plot_loglik,
                                   plot_b0,
                                   plot_b1,
                                   plot_p_x_0,
                                   plot_p_x_1))

print(res)