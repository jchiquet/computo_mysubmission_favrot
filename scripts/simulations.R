

J = 4;   T = 3;   N = 10;   K = 4;

Block = c(1 : K);   Band = c(1 : J);   
Beet = c(1 : N);   DPT = seq(0, 12, length.out = T);     

alpha0 = 0.5;   gamma0 = c(0, - 0.13, - 1.13, - 1.24);       
alpha1 = 0.16;   gamma1 = c(0, 0.24, - 0.14, - 0.15);
sig0 = 1.87;   chi = 0.27;   eta = 0.98      

############################################################################

# Building scenarios #######################################################
# Scenarios "50% Y" and "50% Z" are obtained from scenarios "100% Y" and 
# "100% Z", respectively. For example "50% Y" scenario with 40 trials 
# corresponds to "100% Y" with 20 trials.

seed = 1;   I = 10
data = suppressMessages(simu_data(seed = seed, I = I))

scenarios = list(
  Y = data$dataYZ %>% mutate(Z = NA),
  Z = data$dataYZ %>% mutate(Y = NA),
  YhalfZhalf = data$dataYZ %>% mutate(Y = ifelse(ID <= (I / 2), NA, Y), 
                                      Z = ifelse(ID > (I / 2), NA, Z)),
  YhalfZ = data$dataYZ %>% mutate(Y = ifelse(ID <= (I / 2), NA, Y)),
  YZ = data$dataYZ,
  W = data$dataW
)
#############################################################################


# Inference #################################################################
nadapt = 2000;   niter = 2000

res_inf = NULL

# Inference YZ --------------------------------------------------------------
ID = as.numeric(as.factor(as.character(scenarios$Y$ID)));   
INSEC = as.numeric(as.factor(as.character(scenarios$Y$Insecticide)));   
TIME = as.numeric(scenarios$Y$tscaled);   
ST = as.numeric(as.factor(as.character(scenarios$Y$st)));

K = length(unique(ID));   
L = length(unique(INSEC));   
M = length(unique(ST))

TIME_unique = unique(scenarios$Y$tscaled)[2 : 3]  
T = length(unique(TIME_unique))

for(i in (1 : 5)){
  Y = scenarios[[i]]$Y;   Z = scenarios[[i]]$Z;   
  N = scenarios[[i]]$N;   Q = length(Y)   
  
  data_jags = list(
    "Y" = Y, "Z" = Z, "Q" =  Q, "ID" = ID, "INSEC" = INSEC, 
    "TIME" = TIME, "ST" = ST, "K" = K, "L" = L, "M" = M, 
    "N" = N, "T" = T, "TIME_unique" = TIME_unique
  )
  
  model <- jags.model("Files_for_code/modelYZ.txt", data = data_jags, 
                      n.chains = 2, n.adapt = nadapt)
  
  samples <- coda.samples(model, 
                          variable.names = c("gamma0", "gamma1", "Eff"), 
                          n.iter = niter, thin = 10)
  
  bind = list(samples);   
  names(bind) = paste("samples", names(scenarios)[i], sep = "_")
  
  res_inf = res_inf %>% append(bind)
}

# Inference W ---------------------------------------------------------------
ID = as.numeric(as.factor(as.character(scenarios$W$ID)));   
INSEC = as.numeric(as.factor(as.character(scenarios$W$Insecticide)));   
TIME = as.numeric(scenarios$W$tscaled);   
ST = as.numeric(as.factor(as.character(scenarios$W$st)));
SBIT = as.numeric(as.factor(as.character(scenarios$W$sbit)))

K = length(unique(ID));   L = length(unique(INSEC));   
M = length(unique(ST));   X = length(unique(SBIT))

TIME_unique = unique(scenarios$W$tscaled)[2 : 3]  
T = length(unique(TIME_unique))

W = scenarios$W$W;   Q = length(W)

data_jags = list(
  "W" = W, "Q" =  Q, "ID" = ID, "INSEC" = INSEC, 
  "TIME" = TIME, "ST" = ST,  "SBIT" = SBIT, "K" = K, "L" = L, 
  "M" = M, "X" = X, "T" = T, "TIME_unique" = TIME_unique
)

model <- jags.model("Files_for_code/modelW.txt", data = data_jags, 
                    n.chains = 2, n.adapt = nadapt)

samples <- coda.samples(model, 
                        variable.names = c("gamma0", "gamma1", "Eff"), 
                        n.iter = niter, thin = 10)

res_inf = res_inf %>% append(list("samples_W" = samples))
############################################################################


# Formatting results #######################################################
t = scenarios$Y$tscaled %>% unique
Eff_6_true = (1 - exp(gamma0[2 : J] + gamma1[2 : J] * t[2])) * 100
Eff_12_true = (1 - exp(gamma0[2 : J] + gamma1[2 : J] * t[3])) * 100
truth = c(Eff_6_true, Eff_12_true, gamma0[2 : J], gamma1[2 : J])

n_scenarios = length(scenarios)

esti = lapply(
  res_inf, function(x) summary(x)$statistics %>% as.data.frame %>% 
    rownames_to_column %>% 
    filter(!(grepl("gamma", rowname) & Mean == 0)) %>% 
    select(Mean) %>% as.matrix %>% as.vector
)    

b_inf = lapply(
  res_inf, function(x) summary(x)$quantiles %>% as.data.frame %>% 
    rownames_to_column %>% 
    filter(!(grepl("gamma", rowname) & `2.5%` == 0)) %>% 
    select(`2.5%`) %>% as.matrix %>% as.vector
)    

b_sup = lapply(
  res_inf, function(x) summary(x)$quantiles %>% as.data.frame %>% 
    rownames_to_column %>% 
    filter(!(grepl("gamma", rowname) & `2.5%` == 0)) %>% 
    select(`97.5%`) %>% as.matrix %>% as.vector
)    

parameters = summary(res_inf[[1]])$statistics %>% as.data.frame %>% 
  rownames_to_column %>% 
  filter(!(grepl("gamma", rowname) & Mean == 0)) %>% 
  select(rowname) %>% as.matrix %>% as.vector   

res = do.call(rbind, 
              lapply(c(1 : n_scenarios), 
                     function(x) 
                       data.frame(Truth = truth, Scenario = names(scenarios)[x], 
                                  I = I, seed = seed, 
                                  value = esti[[x]], parameters = parameters, 
                                  b_inf = b_inf[[x]], b_sup = b_sup[[x]]
                       )
              )
)
############################################################################

# save(res, file = paste0("res.I.", I, ".seed.", seed, ".Rdata"))