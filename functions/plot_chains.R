
# function to plot the MCMC algorithm's chains obtained from the inference on real data

plot_chains <- function(samp, nrow = 3){
  
  n = dim(samp[[1]])[1]
  p = dim(samp[[1]])[2]
  
  df_temp = rbind(samp[[1]] %>% as.data.frame %>% mutate(i = c(1 : n), chaine = "1"),
                  samp[[2]] %>% as.data.frame %>% mutate(i = c(1 : n), chaine = "2")) %>% select(- `gamma0[1]`, - `gamma1[1]`)
  
  df_plot = df_temp %>% pivot_longer(cols = c(1 : (p - 2))) %>% as.data.frame %>% 
    mutate(name = recode(name, "alpha0" = "c100", "alpha1" = "c101", "gamma0[2]" = "c102",
                         "gamma0[3]" = "c104", "gamma0[4]" = "c105", "gamma1[2]" = "c106",
                         "gamma1[3]" = "c107", "gamma1[4]" = "c108", "sigma0" = "c109", "chi" = "c110", "eta" = "c111",
                         "Eff[2,1]" = "c112", "Eff[3,1]" = "c113", "Eff[4,1]" = "c114", 
                         "Eff[2,2]" = "c115", "Eff[3,2]" = "c116", "Eff[4,2]" = "c117"))
  
  titre1 = TeX("$\\alpha_0$");   titre2 = TeX("$\\alpha_1$");   
  
  titre3 = TeX("$\\gamma_0 - Mavrik Jet$");   titre4 = TeX("$\\gamma_0 - Movento$");   titre5 = TeX("$\\gamma_0 - Teppeki$");   
  titre6 = TeX("$\\gamma_1 - Mavrik Jet$");   titre7 = TeX("$\\gamma_1 - Movento$");   titre8 = TeX("$\\gamma_1 - Teppeki$");   
  
  titre9 = TeX("$\\sigma_0$");   titre10 = TeX("$\\chi$");   titre11 = TeX("$\\eta$");   
  
  titre12 = TeX("$Ef_6 - Mavrik Jet$");   titre13 = TeX("$Ef_6 - Movento$");    titre14 = TeX("$Ef_6 - Teppeki$")
  titre15 = TeX("$Ef_{12} - Mavrik Jet$");   titre16 = TeX("$Ef_{12} - Movento$");    titre17 = TeX("$Ef_{12} - Teppeki$")
  
  df_plot = df_plot %>% mutate(name = as.factor(name))
  
  levels(df_plot$name) = c(titre1, titre2, titre3, titre4, titre5, titre6, titre7, titre8, titre9, titre10, titre11,
                           titre12, titre13, titre14, titre15, titre16, titre17)
  
  ggplot(df_plot) + geom_line(aes(x = i, y = value, color = chaine)) + facet_wrap(~ name, scales = "free_y", nrow = nrow, labeller = label_parsed) + theme(legend.position = "none") + xlab("Iteration") + ylab("Sampled value")
}