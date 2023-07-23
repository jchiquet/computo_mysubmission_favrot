simu_data <- function(seed, I){
 
 set.seed(seed)
 
 ID = c(1 : I);   
 
 data = expand_grid(ID = c(1 : I),   Block = c(1 : K),   Band = c(1 : J), 
                    Beet = c(1 : N),   DPT = seq(0, 12, length.out = T))
 
 data = data %>% 
  mutate(Insecticide = ifelse(DPT <= 0, "T 1", paste("T", Band)),
         tscaled = scale(data$DPT)) %>% 
  mutate(st = paste(ID, Insecticide), 
         sbit = paste(ID, Block, Insecticide, DPT))
 
 Insecticide = data$Insecticide %>% unique %>% sort
 names(gamma0) = Insecticide
 names(gamma1) = Insecticide
 
 beta0 = rnorm(I, sd = sig0);   
 u = rnorm(I * J, sd = chi);   
 epsi = rnorm(I * J * K * (T - 1) + I * K, sd = eta);   
 
 names(beta0) = data$ID %>% unique;   
 names(u) = data$st %>% unique
 names(epsi) = data$sbit %>% unique
 
 data = data %>% mutate(alpha0 = alpha0,
                        alpha1 = alpha1,
                        N = N,
                        beta0 = recode(ID, !!!beta0),
                        gamma0 = recode(Insecticide, !!!gamma0), 
                        gamma1 = recode(Insecticide, !!!gamma1),
                        u = recode(st, !!!u),
                        epsi = recode(sbit, !!!epsi))
 
 data$lb = exp(data$alpha0 + data$beta0 + data$gamma0 + (data$alpha1 + 
                                                          data$gamma1) * data$tscaled + data$u + data$epsi)
 
 data$W = sapply(c(1 : (I * J * K * T * N)), 
                 function(x) rpois(1, data$lb[x]))
 
 dataYZ = data %>% group_by(ID, Block, Band, DPT) %>% 
  summarise(Insecticide = unique(Insecticide), 
            tscaled = unique(tscaled), st = unique(st), N = unique(N), 
            Y = sum(W), Z = sum(W > 0)) %>% as.data.frame 
 
 dataW = data %>% select(- alpha0, - alpha1, - beta0, - gamma0, - gamma1, 
                         - u, - epsi, - lb)
 
 return(list("dataYZ" = dataYZ, "dataW" = dataW))
}                        
