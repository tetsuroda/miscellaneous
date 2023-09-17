
# ###########################
# Basic functions
# ###########################

logit <- function(p){
  return(
    log(p/(1-p))
  )
}

fun_beta0 <- 
  function(x1,x2,beta1, beta2) {
    eta_wo_intercept <- beta1*x1 + beta2*x2
    
    # adjust beta0 to get the assumed prevalence
    beta0 <- logit(prevalence) - eta_wo_intercept 
    return(beta0)
  }


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7949469/
calc_ppvs <-
  function(.prevalence, .sens, .spe) {
    
    ppv <- .prevalence*.sens/(.prevalence*.sens + (1-.prevalence)*(1-.spe))
    
    return(ppv) 
  }


# ###########################
# Create-data function 
# ###########################
# Assume we have a binary outcome Y with two binary covariates, x1 and x2, in the population model.
# Using pre-specified beta1 and 2 as well as the proportion of covariates,
# the function produces a toy dataset.

make_data <- function(.prevalence,
                      .beta1,
                      .beta2,
                      .N,
                      .x1_prevalence,
                      .x2_prevalence) {
  
  prevalence <- .prevalence
  beta1 <- .beta1
  beta2 <- .beta2

  N <- .N
  x1_prevalence <- .x1_prevalence
  x2_prevalence <- .x2_prevalence
  
  x1 <- as.numeric(runif(n=N) <= x1_prevalence)
  x2 <- as.numeric(runif(n=N) <= x2_prevalence)
  
  # Get a population intercept parameter
  # as a mean of individual intercepts 
  beta0 <-
    data.frame(x1,x2,beta1 = beta1, beta2 = beta2) %>%
    pmap(fun_beta0) %>%
    unlist() %>% mean()
  
  # Calculate individual probability of getting an event
  individual_p <-
    plogis(beta1*x1 + beta2*x2+beta0)
  
  # Check the mean of individual probability similar to the pre-specified prevalence
  # mean(individual_p)
  
  df <-
    data.frame(x1,x2,individual_p)
  
  # Return the dataset with the population beta0
  return(list(dataset=df, beta0=beta0))
  
}

# ###########################
# Simulation function
# ###########################
# The function produces the random FN or FP in the outcome of the original sample
# and run logistic regression in the original sample and the modified sample
# to see if misclassification in the outcome affects the OR estimates
# NB: Inside the function, {caret} is used.

sim_results <-
  function(.reps, .data, .sens, .spe, .prevalence) {
    
    # parameters
    sens <- .sens 
    spe <- .spe　
    prevalence <- .prevalence
    N <- nrow(.data)
    
    # true y
    .data$y <- rbinom(n=N, size = 1, prob = .data$individual_p)
    
    # can confirm the estimated prevalence
    # empirical_prevalence <- mean(.data$y)
    
    # Theoretical ppv
    ppv <- round(calc_ppvs(prevalence, sens, spe),3)
    
    .data$y_def1[.data$y==1] <- rbinom(n=sum(.data$y==1), size = 1, prob = sens)
    .data$y_def1[.data$y==0] <- (rbinom(n=sum(.data$y==0), size = 1, prob = spe)-1)*-1
    
    # Logistic regression
    fit1 <- glm(y~x1+x2, family = "binomial", data = .data)
    fit2 <- glm(y_def1~x1+x2, family = "binomial", data = .data)
    
    # 2by2 table with the y_def1
    out_table <- table(.data$y_def1,.data$y)[c(2,1),c(2,1)]
    
    # Output the results in a tidy format
    res <-
      data.frame(
        model= "trueY",
        broom::tidy(fit1,  exponentiate = FALSE)[,1:3]
      ) %>% bind_rows(
        data.frame(
          model= paste0("sens_",sens,"_spe_",spe,"_ppv_", ppv),
          broom::tidy(fit2, exponentiate = FALSE)[,1:3],
          insample_sens = caret::sensitivity(out_table),
          insample_spe = caret::specificity(out_table),
          insample_ppv = caret::posPredValue(out_table),
          TP = out_table[1,1],
          FN = out_table[2,1],
          FP = out_table[1,2],
          TN = out_table[2,2]
        ) 
      )
    
    return(res)
  }

# ###########################
# Visualize function
# ###########################
make_plots <-
  function(.input,
           .prevalence,
           .x1_prevalence,
           .x2_prevalence,
           .reps) {
    
    
    prevalence <- .prevalence
    x1_prevalence <- .x1_prevalence
    x2_prevalence <- .x2_prevalence
    reps <- .reps
    
    df_plot <-
      .input %>%
      group_by(model, term) %>%
      mutate(beta_low25 = quantile(estimate, probs=0.025),
             beta_mean = median(estimate),
             beta_high975 = quantile(estimate, probs=0.975),
             bias_betas = mean(diff_betas)
             ) %>%
      distinct(beta_low25,beta_mean, beta_high975, true_betas, bias_betas) %>%
      mutate(labels = paste(term, model)) %>%
      ungroup()
    
    res_plot <-
      df_plot %>%
      group_by(model, term) %>%
      ggplot(aes(y=labels)) +
      geom_point(aes(x=exp(beta_mean))) +
      geom_errorbarh(height=0.25,aes(xmin=exp(beta_low25), xmax=exp(beta_high975))) +
      geom_point(aes(x = exp(true_betas), y = labels), col = "red", shape = 2) +
      theme_bw() +
      labs(x="exp(betas)[odds]") +
      ggtitle(paste0("Pr(Y=1)=",prevalence,", N=",nrow(df),", reps=",reps),
              subtitle = paste0("Pr(X1=1)=",x1_prevalence,", Pr(X2=1)=",x2_prevalence))
    
    res_plot2 <-
      .input %>%
      mutate(labels = paste(term, model)) %>%
      ggplot() +
      geom_jitter(aes(y=labels, x=exp(estimate)),alpha = 0.3) +
      geom_boxplot(aes(y=labels,x=exp(estimate)), alpha = 0, col="blue") +
      geom_point(data = .input %>% distinct(term, model, true_betas) %>% mutate(labels = paste(term, model)),
                 aes(x = exp(true_betas), y = labels), col = "red", shape = 2) +
      theme_bw() +
      labs(x="exp(betas)[odds]") +
      ggtitle(paste0("Pr(Y=1)=",prevalence,", N=",nrow(df),", reps=",reps),
              subtitle = paste0("Pr(X1=1)=",x1_prevalence,", Pr(X2=1)=",x2_prevalence))
    
    
    return(list(data_for_plot=df_plot, 
                errorbarplot=res_plot,
                boxplot=res_plot2))
  }


# ###########################
# Visualize PPV
# ###########################

df_xy <-
  expand.grid(seq(0,1,0.01),
              seq(0,1,0.01))
colnames(df_xy) <- c(".sens",".spe")

df_ppvs <-
  data.frame(.prevalence = prevalence,
             df_xy) %>%
  mutate(ppv = calc_ppvs(.prevalence, .sens, .spe))

df_ppvs %>%
  ggplot(aes(x=.sens, y=.spe, fill=ppv)) +
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  ggtitle(paste0("Pr(Y=1)=",prevalence))


# ###########################
# Execute
# ###########################

# Create a dataset
prevalence <- 0.2
odds1 <- 3; beta1 <- log(odds1)
odds2 <- 3; beta2 <- log(odds2)
N <- 1000
x1_prevalence <- 0.3
x2_prevalence <- 0.8

df_out <- make_data(.prevalence = prevalence,
                .beta1 = beta1,
                .beta2 = beta2,
                .N = N,
                .x1_prevalence = x1_prevalence,
                .x2_prevalence = x2_prevalence)

df <- df_out$dataset
beta0 <- df_out$beta0

# NB: We could repeat dataset sampling as well, but
# I wanted to see only the effect of misclassification here.

# Scenario1
# #################################################
# As long as Spe = 1, random misclassification does not reduce ensitivity in RR.
# However there should be some bias in OR

reps <- 1000
sens <- 0.8
spe <- 0.99 #If 1 no FP is produced and no stochastic results to be out

# Run logistic regression
# Expect time needed as follows with reps = 1000
# user  system elapsed 
# 16.471   0.092  16.577 
system.time(
out <-
  seq(1:reps) %>%
  map(sim_results, .data = df, .sens = sens, .spe = spe, .prevalence = prevalence) %>%
  bind_rows() %>%
  mutate(true_betas = case_when(
    term == "x1" ~ beta1,
    term == "x2" ~ beta2,
    TRUE ~ beta0
  )) %>%
  mutate(diff_betas = estimate - true_betas)
)

# Make plots
out_plot <-
  make_plots(.input=out,
          .prevalence = prevalence,
          .x1_prevalence = x1_prevalence,
          .x2_prevalence = x2_prevalence,
          .reps = reps)

out_plot$errorbarplot
out_plot$boxplot

# Calculate bias in percent
out_plot$data_for_plot %>%
  mutate(bias_betas_percent = round(bias_betas/true_betas*100,2)) %>%
  select(term, model, bias_betas_percent)
  

# Scenario 2
# #################################################
# Unless spe = 1, even sensitivity = 1 produces bias

reps <- 1000
sens <- 0.99
spe <-0.9

out2 <-
  seq(1:reps) %>%
  map(sim_results, .data = df, .sens = sens, .spe = spe, .prevalence = prevalence) %>%
  bind_rows() %>%
  mutate(true_betas = case_when(
    term == "x1" ~ beta1,
    term == "x2" ~ beta2,
    TRUE ~ beta0
  )) %>%
  mutate(diff_betas = estimate - true_betas)

out_plot2 <-
  make_plot(.input=out2,
            .prevalence = prevalence,
            .x1_prevalence = x1_prevalence,
            .x2_prevalence = x2_prevalence,
            .reps = reps)

out_plot2$errorbarplot
out_plot2$boxplot

out_plot2$data_for_plot %>% 
  mutate(bias_betas_percent = round(bias_betas/true_betas*100,2)) %>%
  select(term, model, bias_betas_percent)


  
# Scenario 3
# #################################################
# See what happens with a realistic setting

reps <- 1000
sens <- 0.8
spe <-0.8

out3 <- 
  seq(1:reps) %>%
  map(sim_results, .data = df, .sens = sens, .spe = spe, .prevalence = prevalence) %>%
  bind_rows() %>%
  mutate(true_betas = case_when(
    term == "x1" ~ beta1,
    term == "x2" ~ beta2,
    TRUE ~ beta0
  )) %>%
  mutate(diff_betas = estimate - true_betas)

out_plot3 <-
  make_plot(.input=out3,
            .prevalence = prevalence,
            .x1_prevalence = x1_prevalence,
            .x2_prevalence = x2_prevalence,
            .reps = reps)

out_plot3$errorbarplot
out_plot3$boxplot
  
out_plot3$data_for_plot %>% 
  mutate(bias_betas_percent = round(bias_betas/true_betas*100,2)) %>%
    select(term, model, bias_betas_percent)
  
  