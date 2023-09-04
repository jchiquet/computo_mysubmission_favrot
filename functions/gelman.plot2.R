
source(file = "functions/gelman.R")

# function adapted from gelman.plot, using ggplot instead of plot

gelman.plot2 <- function (x, coef, bin.width = 10, max.bins = 50, confidence = 0.95, 
                          transform = FALSE, autoburnin = TRUE, auto.layout = TRUE, 
                          ask, col = 1:2, lty = 1:2, xlab = "last iteration in chain", 
                          ylab = "shrink factor", type = "l", ncol = 4, ...){
  exclude = which((summary(x)$statistics %>% rownames) %in% c("gamma0[1]", "gamma1[1]"))
  x = x[, - exclude]
  if (missing(ask)) {
    ask <- if (is.R()) {
      dev.interactive()
    }
    else {
      interactive()
    }
  }
  x <- as.mcmc.list(x)
  oldpar <- NULL
  on.exit(par(oldpar))
  
  y <- gelman.preplot(x, bin.width = bin.width, max.bins = max.bins, 
                      confidence = confidence, transform = transform, autoburnin = autoburnin)
  all.na <- apply(is.na(y$shrink[, , 1, drop = FALSE]), 2, 
                  all)
  l = list()
  df_value = data.frame(median = NULL, bsup = NULL, Param = NULL, iter = NULL)
  coef = rownames(summary(x)$statistics)
  if (!any(all.na)){ 
    for (j in 1:nvar(x)) {
      df_temp = as.data.frame(y$shrink[, j, ]);   df_temp$Param = coef[j];    df_temp$iter = y$last.iter;   colnames(df_temp)[c(1, 2)] = c("median", "97.5%")
      df_value = rbind(df_value, df_temp)
    }
  }    
  
  df_res = df_value %>% pivot_longer(cols = c(median, `97.5%`)) %>% as.data.frame %>% 
    mutate(Param = recode(Param, "alpha0" = "c100", "alpha1" = "c101", "gamma0[2]" = "c102",
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
  
  df_res = df_res %>% mutate(Param = as.factor(Param))
  
  levels(df_res$Param) = c(titre1, titre2, titre3, titre4, titre5, titre6, titre7, titre8, titre9, titre10, titre11,
                           titre12, titre13, titre14, titre15, titre16, titre17)
  
  return(ggplot(df_res) + geom_line(aes(x = iter, y = value, color = name, linetype = name)) + facet_wrap(~ Param, ncol = ncol, labeller = label_parsed) +
           xlab("Last iteration in chain") + ylab("Shrink factor") + theme(legend.position = "bottom") + theme(legend.title = element_blank()) +
           scale_linetype_manual(values = c("dotted", "solid")) + 
           scale_color_manual(values = c("red", "black")))
}