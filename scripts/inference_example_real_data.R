

# Jags code for the model  #################################################
modelstringYZ = "
    model {

        # Likelihood #####################################################
        for (i in 1:Q){
            Y[i] ~ dpois(N[i] * lb[i])
            Z[i] ~ dbinom(pi[i], N[i])

            log(lb[i]) = beta0[ID[i]] + gamma0[INSEC[i]] + (alpha1 + 
                         gamma1[INSEC[i]]) * TIME[i] + u[ST[i]] + epsi[i]

            pi[i] = 1 - exp(- lb[i])
            epsi[i] ~ dnorm(0, pi_eps)
        }

        for (j in 1:K){
            beta0[j] ~ dnorm(alpha0, tau0)
        }

        for (c in 1:M){
            u[c] ~ dnorm(0, invchi)
        }

        gamma0[1] = 0
        gamma1[1] = 0

        # Priors #########################################################
        for (s in 2:L){
            gamma0[s] ~ dnorm(0, 0.001)
            gamma1[s] ~ dnorm(0, 0.001)
        }

        alpha0 ~ dnorm(0, 0.001)
        alpha1 ~ dnorm(0, 0.001)
        sigma0 ~ dunif(0, 10)
        chi ~ dunif(0, 10)
        eta ~ dunif(0, 10)

        # Derived Quantities #############################################
        tau0 = pow(sigma0, -2)
        invchi = pow(chi, -2)
        pi_eps = pow(eta, -2)

        for (h in 2:L){
            for(t in 1 : T){
                Eff[h, t]  = (1 - exp(gamma0[h] + gamma1[h] * 
                             TIME_unique[t])) * 100
            }
        }  
    }
"

writeLines(modelstringYZ, con  =  "Files_for_code/modelYZ.txt")
############################################################################


# Inference example on the extract of the real dataset #####################
load(file = "Files_for_code/real_data_extract.Rdata")
data = real_data_extract %>% 
  mutate(tscaled = scale(DPT), st = paste(ID, Insecticide))

# Building scenarios -------------------------------------------------------
scenarioY = data %>% mutate(Z = NA)
scenarioYhalfZhalf = data %>% mutate(Y = ifelse(ID == "2020 - B1A97", Y, NA), 
                                     Z = ifelse(ID == "2020 - B1A97", NA, Z))
scenarioYhalf = data %>% filter(ID == "2020 - B1A97") %>% mutate(Z = NA)
scenarioZhalf = data %>% filter(ID == "2020 - u1CwE") %>% mutate(Y = NA)
# ---------------------------------------------------------------------------

data = scenarioYhalfZhalf

Y = data$Y;   Q = length(Y);   N = data$N;   Z = data$Z

ID = as.numeric(as.factor(as.character(data$ID)));   
INSEC = as.numeric(as.factor(as.character(data$Insecticide)));   
TIME = as.numeric(data$tscaled);   
ST = as.numeric(as.factor(as.character(data$st)));   

K = length(unique(ID));   L = length(unique(INSEC));   
M = length(unique(ST));   

df_TIME = suppressMessages(data %>% 
                             group_by(DPT, tscaled) %>% 
                             summarise(n = n()) %>% as.data.frame)

TIME_unique = approx(df_TIME$DPT, df_TIME$tscaled, xout = c(6, 12))$y;   
T = length(unique(TIME_unique))

data_jags = list(
  "Y" = Y, "Z" = Z, "Q" =  Q, "ID" = ID, "INSEC" = INSEC, 
  "TIME" = TIME, "ST" = ST, "K" = K, "L" = L, "M" = M, 
  "N" = N, "T" = T, "TIME_unique" = TIME_unique
)

nadapt = 2000;   niter = 2000

model <- jags.model("Files_for_code/modelYZ.txt", data = data_jags, 
                    n.chains = 2, n.adapt = nadapt)

samples <- coda.samples(model, 
                        variable.names = c("gamma0", "gamma1", "Eff"), 
                        n.iter = niter, thin = 10)
############################################################################