# ##############################################################################
# Trade-offs between accuracy measures for electronic healthcare data algorithms
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3264740/#__ffn_sectitle
# #############################################################################

# visualization function for bias %
# ################################################
make_bias_percent <-
  function(sens, spe, .prop_exposed, .I_e, .I_u) {
    prop_exposed <- .prop_exposed
    I_e <- .I_e
    I_u <- .I_u
    RR_true <- round(I_e/I_u,3)
    a <- sens*prop_exposed*I_e + (1-spe)*prop_exposed*(1-I_e)
    b <- (1-sens)*prop_exposed*I_e + spe*prop_exposed*(1-I_e)
    c <- sens*(1-prop_exposed)*I_u + (1-spe)*(1-prop_exposed)*(1-I_u)
    d <- (1-sens)*(1-prop_exposed)*I_u + spe*(1-prop_exposed)*(1-I_u)
    
    a_true <- prop_exposed*I_e
    b_true <- prop_exposed*(1-I_e)
    c_true <- (1-prop_exposed)*I_u
    d_true <- (1-prop_exposed)*(1-I_u)
    
    RR_obs <- (a/(a+b))/(c/(c+d))
    bias_multiplier <- round(RR_obs/RR_true,3)
    bias_percent <- round(100*(RR_obs - RR_true)/RR_true,3)
    
    out <- data.frame(a,b,c,d,
                      a_true,b_true,c_true,d_true,
                      RR_obs, RR_true,
                      bias_multiplier,bias_percent)
    return(out)
    
  }


# Execute
# ##################################################
df_se_sp <-
  expand.grid(seq(0,1,0.01),
            seq(0,1,0.01))
colnames(df_se_sp) <- c("sens","spe")

# fixed parameters
prop_exposed <- 0.7
I_u <- 0.1
I_e <- 0.2
(incidence_y <- prop_exposed*I_e + (1-prop_exposed)*I_u)
(RR_true <- round(I_e/I_u,3))

chubak_out <-
  df_se_sp %>%
  pmap_df(make_bias_percent,
       .prop_exposed = prop_exposed,
       .I_u = I_u,
       .I_e = I_e) 

N <- 1000

chubak_out2 <-
  df_se_sp %>% bind_cols(chubak_out) %>%
  mutate(a = as.integer(a*N),
         b = as.integer(b*N),
         c = as.integer(c*N),
         d = as.integer(d*N),
         a_true = as.integer(a_true*N),
         b_true = as.integer(b_true*N),
         c_true = as.integer(c_true*N),
         d_true = as.integer(d_true*N),
  )

chubak_out2 %>%
  ggplot(aes(x=sens, y=spe, fill=bias_percent)) +
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  ggtitle(paste0("Pr(Y=1)=",incidence_y,", Pr(exposed)=",prop_exposed,", true RR=",RR_true),
          subtitle=paste0("Incidence in exposeds=",I_e,", Incidence in unexposeds=",I_u))
