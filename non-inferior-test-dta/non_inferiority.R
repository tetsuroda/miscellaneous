# ##################################################
# Non-inferiority test for sensitivity of diagnosis
# #################################################

tab <- matrix(c(555,45,40,860), ncol = 2, byrow = T)
marg <- -0.01 # margin under H0 
n <- sum(tab)
alpha = 0.05

x11_obs = tab[1,1]
x10_obs = tab[1,2] # new + & ref -
x01_obs = tab[2,1] # new - & ref + 
x00_obs = tab[2,2]

p01_obs = x01_obs / n
p10_obs = x10_obs / n

sen1_obs = (x11_obs + x10_obs) / n # new diagnosis
sen2_obs = (x11_obs + x01_obs) / n # reference
theta_obs =  sen1_obs - sen2_obs

# Nam (1997) https://www.jstor.org/stable/2533508?origin=JSTOR-pdf
# -----------------------------------------------------------------
a <- 2*n
b <- (2*n + x01_obs - x10_obs)*marg - (x01_obs + x10_obs)
c <- -x01_obs*marg*(1-marg)
p01_mle <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
q01_mle <- 1- p01_mle
p10_mle <- p01_mle + marg

z <- ((x10_obs*q01_mle - (n - x01_obs)*p10_mle) / (p10_mle*(1 - p01_mle - p10_mle))) * sqrt((p10_mle + p01_mle - marg^2)/n)
# p-value
pvalue_nam <- pnorm(z, lower.tail = F)

# Lower limit of CI
CI_l_nam <- theta_obs - qnorm(1-alpha)*sqrt((p10_mle + p01_mle - marg^2)/n)

out_nam <- 
data.frame(type = 'Nam',
           n_matched_pairs = n, 
           new_sens = sen1_obs, 
           reference_sens = sen2_obs,
           diff_new_from_ref = theta_obs,
           lower_90_CI = CI_l_nam,
           inferioirty_p_value = pvalue_nam)


# Liu et al (2002) https://pubmed.ncbi.nlm.nih.gov/11782062/
# -----------------------------------------------------------------
c = x10_obs - x01_obs
d = x01_obs + x10_obs

# If this is TRUE, the hypothesis is rejected = non-inferiority is established
# also according to Liu, 1-alpha/2 is recommended in ICH E9
z_l = (c+n*abs(marg))/sqrt(d-n*theta_obs^2)
z_l >= qnorm(1-alpha/2) 

# p-value
pvalue_liu = pnorm(z_l, lower.tail = F)

# If this TRUE, the hypothesis is rejected = non-superiority is established
z_u = (c-n*abs(marg))/sqrt(d-n*theta_obs^2)
z_u <= -1*qnorm(1-alpha/2) 

# confidence intervals
# sigma_obs is defined as in Liu
# significance level is explained in Walker_Nowacki_2010
sigma_obs = sqrt(((p01_obs+p10_obs)-theta_obs^2)/n)
(CI_l_liu = theta_obs - qnorm(1-alpha)*sigma_obs) 

# when to test equivalence by two one sided test (TOST), 
# (1-2*alpha)100% CI is used
# (CI_u_liu = theta_obs + qnorm(1-2*alpha)*sigma_obs) 


out_liu <-
data.frame(type = 'Liu',
           n_matched_pairs = n, 
           new_sens = sen1_obs, 
           reference_sens = sen2_obs,
           diff_new_from_ref = theta_obs,
           lower_95_CI = CI_l_liu,
           inferioirty_p_value = pvalue_liu)


# Results
# -----------------------------------------------------------------
out_nam %>%
  bind_rows(out_liu)




