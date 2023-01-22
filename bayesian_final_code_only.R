#----------------------------
# Dickson, 6369693, only code
#----------------------------

options(scipen=999)
library(tidyverse)
library(ggpubr)
library(sn)
library(kableExtra)
library(bain)

sleep <- read.csv("sleephealth.csv")                             # import the data
sleep <- sleep[,c(6,7,8,22)]                                     # keep these variables
names(sleep) <- c("active", "duration", "onset", "screen")       # rename the columns
sleep$interaction <- sleep$onset*sleep$screen
sleep <- as.data.frame(scale(sleep, center = F))                 # scale/center the data 
y  <- sleep$duration     # assign outcome 
x1 <- sleep$active    # assign predictors
x2 <- sleep$screen
x3 <- sleep$onset
x4 <- sleep$interaction
n  <- nrow(sleep)

# initial values
beta0 <- 0
beta1 <- 0.07
beta2 <- 0.06
beta3 <- 0.05
sigma2 <- 0.5
step <- 0.0005

Gibbs_MH_Sampler <- function(y, x1, x2=NULL, x3=NULL, n, mu.00, sigma2.00, mu.10, sigma2.10, mu.30, sigma2.30, A.0, B.0){
  # Iterations----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  burnIN <- 30000
  thin <- 1000
  nITER <- burnIN + thin
  # Storage for the estimates-------------------------------------------------------------------------------------------------------------------------------------------------------
  beta0s <- c()
  beta1s <- c()
  beta2s <- c()
  beta3s <- c()
  sigma2s<- c()
  js<- c()
  # Start sampling----------------------------------------------------------------------------------------------------------------------------------
  for (chain in 1:2){
    j <- 0 
    # Iterate over parameters-------------------------------------------------------------------------------------------------------------------------
    for(iter in 1:nITER){
      j <- j+1 
      # Sample beta0
      beta0 <- rnorm(1, mean = ((sum(y - beta1*x1 -beta2*x2 - beta3*x3)/sigma2) + mu.00/sigma2.00) / (n/sigma2 + 1/sigma2.00), sd = sqrt(1/(n/sigma2 + 1/sigma2.00)))
      
      # Sample beta1
      beta1 <- rnorm(1, mean = ((sum(x1*(y - beta0 - beta2*x2 - beta3*x3))/sigma2) + mu.10/sigma2.10) / (sum(x1^2)/sigma2 + 1/sigma2.10), sd = sqrt(1/(n/sigma2 + 1/sigma2.00)))
      
      # Sample beta2 
      #prop  <- rnorm(1, mean = beta2, sd = step) 
      prop  <- rsn(xi=1.256269, omega=1.605681, alpha=5) # skew-normal distribution
      ref  <- runif(1)     
      poa <- exp((dt(prop, df=2, log = T))) / exp((dt(beta2, df=2, log = T))) * (dnorm(beta2, mean = beta2, sd= step) / dnorm(prop, mean=beta2, sd=step))
      accept  <- ifelse(poa >= ref, 1, 0)                                         
      if (accept == 1){
        beta2 <- prop}
      
      # Sample beta_3
      beta3 <- rnorm(1, mean = ((sum(x3*(y - beta0 - beta2*x2 - beta1*x1))/sigma2) + mu.30/sigma2.30) / (sum(x3^2)/sigma2 + 1/sigma2.30), sd = 1/ (sum(x3^2)/sigma2 + 1/sigma2.30))
      
      # Sample sigma^2 
      yhat <-  beta0 + beta1*x1 + beta2*x2 + beta3*x3
      A <- A.0 + n/2
      B <- sum((y-yhat)^2) /2 + B.0
      sigma2   <-  1 /rgamma(1, A, B)
      
      # Storage vector for estimates--------------------------------------------------------------------------------------------------------------------
      #keep[i,] <- c(beta0, beta1, beta2, beta3, sigma2)
      beta0s[iter]  <- beta0 
      beta1s[iter]  <- beta1
      beta2s[iter]  <- beta2
      beta3s[iter]  <- beta3
      sigma2s[iter]  <- sigma2
      js[iter] <- j
    }
    
    # Store the results for each chain----------------------------------------------------------------------------------------------------------------
    if (chain == 1){                                                                 # CHAIN 1      
      dat <- cbind(beta0s, beta1s, beta2s, beta3s, sigma2s, chain, js)
      chain1 <- as.data.frame(dat[(burnIN+1):nITER,])}
    else{                                                                           # CHAIN 2
      dat <- cbind(beta0s, beta1s, beta2s, beta3s, sigma2s, chain, js)
      chain2 <- as.data.frame(dat[(burnIN+1):nITER,])
      pooled_chains <- rbind(chain1, chain2)                                          # POOLED CHAINS 
    }
  }
  # Sample statistics for each chain----------------------------------------------------------------------------------------------------------------
  chain1_sample <- as.data.frame(matrix(0, 5, 5))                                   # CHAIN 1                
  rownames(chain1_sample)<- c("Intercept", "beta1", "beta2", "beta3", "sigma")
  colnames(chain1_sample)<- c("Mean", "SD", "CI 2.5%", "CI 97.5", "MC Error")
  chain1b <- chain1 %>% select(-js, -chain)
  chain1_sample[,1]<- apply(chain1b, 2, mean)
  chain1_sample[,2]<- apply(chain1b, 2, sd)
  chain1_sample[,3]<- apply(chain1b, 2, quantile, 0.025)
  chain1_sample[,4]<- apply(chain1b, 2, quantile, 0.975)
  sds <- chain1_sample[,2]
  chain1_sample[,5]<- sds[1:5]/sqrt(nITER)
  
  chain2_sample <- as.data.frame(matrix(0, 5, 6))                                   # CHAIN 2             
  rownames(chain2_sample)<- c("Intercept", "beta1", "beta2", "beta3", "sigma")
  colnames(chain2_sample)<- c("Mean", "SD", "CI 2.5%", "CI 97.5", "MC Error")
  chain2b <- chain2 %>% select(-js, -chain)
  chain2_sample[,1]<- apply(chain2b, 2, mean)
  chain2_sample[,2]<- apply(chain2b, 2, sd)
  chain2_sample[,3]<- apply(chain2b, 2, quantile, 0.025)
  chain2_sample[,4]<- apply(chain2b, 2, quantile, 0.975)
  sds <- chain2_sample[,2]
  chain2_sample[,5]<- sds[1:5]/sqrt(nITER)
  
  pooled_sample <- as.data.frame(matrix(0, 5, 6))                                   # POOLED CHAINS
  rownames(pooled_sample)<- c("beta0", "beta1", "beta2", "beta3", "sigma")
  colnames(pooled_sample)<- c("Mean", "SD", "CI 2.5%", "CI 97.5", "MC Error")
  pooled_chainsb <- pooled_chains %>% select(-js,-chain) 
  pooled_sample[,1]<- apply(pooled_chainsb, 2, mean)
  pooled_sample[,2]<- apply(pooled_chainsb, 2, sd)
  pooled_sample[,3]<- apply(pooled_chainsb, 2, quantile, 0.025)
  pooled_sample[,4]<- apply(pooled_chainsb, 2, quantile, 0.975)
  sds <- pooled_sample[,2]
  pooled_sample[,5]<- sds[1:5]/sqrt(nITER)
  
  # DIC---------------------------------------------------------------------------------------------------------------------------------------------
  
  two_ll <- rep(NA, nrow(pooled_sample))
  for(i in 1:nrow(pooled_sample)){
    two_ll[i] <- -2*sum(dnorm(y, mean = pooled_chains[i,2] + pooled_chains[i,3]*x1 + pooled_chains[i,4]*x2 + pooled_chains[i,5]*x3, sd =sqrt(pooled_chains[i,6]),log=T))}
  d_bar = mean(two_ll)
  d_hat = 2*sum(dnorm(y, mean = pooled_sample$Mean[1] + pooled_sample$Mean[2]*x1 + pooled_sample$Mean[3]*x2 + pooled_sample$Mean[4]*x3, sd = sqrt(pooled_sample$Mean[5]), log=T))
  
  DIC <- d_hat + 2*(d_bar-d_hat)
  
  # Return these objects----------------------------------------------------------------------------------------------------------------------------
  my_list1 <- list("chain1"=chain1, "chain1_sample"=chain1_sample, "chain2"=chain2, "chain2_sample"=chain2_sample, "pooled_chains"=pooled_chains, "pooled_sample"=pooled_sample,"DIC"=DIC, "nITER"=nITER, "thin"=thin, "burnIN"=burnIN)
  list2env(my_list1 ,.GlobalEnv) # to the environment only 
  my_list2 <- list(chain1_sample, chain2_sample, pooled_sample)
  return(my_list2)
}

# Model 1
mod_1 <- Gibbs_MH_Sampler(y=sleep$duration,x1=sleep$active,x2=sleep$screen,x3=sleep$interaction,n=168,mu.00=0,sigma2.00=0.25,mu.10=0,sigma2.10=0.5,mu.30=0,sigma2.30=0.5,A.0=0.00001,B.0=0.00001)

# Posterior predictive check
simdata <- matrix(0,nrow=nrow(pooled_chains),ncol=3)   # storage for simulated data
obsdata <- matrix(0,nrow=nrow(pooled_chains), ncol=3)  # storage for observed data
vars <- cbind(rep(1,n), x1,x2,x3)                      # observed data matrix 

# Predicted data, test-statistic (outliers), discrepancy measure (skew)
for (i in 1:nrow(pooled_chains)){
  betas<- as.numeric(pooled_chains[i, 2:5])           # pooled betas
  sd<- sqrt(pooled_chains[i,6])                       # pooled standard deviation
  R_est<- vars%*%betas                                # predicted y (sleep minutes)
  R_fit<- R_est + rnorm(n,0,sd)                       # dependent variable
  simres<- R_fit - R_est                              # residual errors
  simdata[i,1]<- max(R_est)-min(R_est)                                                        # Outliers (range)
  simdata[i,2]<- (sum((simres-mean(simres))^3)/n)/(sum((simres-mean(simres))^2)/n)^(3/2)      # Skew
  simdata[i,3]<- cor(simres, R_fit)                                                           # crude measure of homoscedasticity
}

# Predicted data, test-statistic (outliers), discrepancy measure (skew)
for (i in 1:nrow(pooled_chains)){
  betas<- as.numeric(pooled_chains[i, 2:5])     
  sd<- sqrt(pooled_chains[i,6])
  R_est<- vars%*%betas
  R_obs<- y
  obsres <- R_obs - R_est
  obsdata[i,1]<- max(R_obs)-min(R_obs)
  obsdata[i,2]<- (sum((obsres-mean(obsres))^3)/n)/(sum((obsres-mean(obsres))^2)/n)^(3/2)
  obsdata[i,3]<- cor(obsres, R_obs)                                                           # crude measure of homoscedasticity
}

# Posterior predictive p-value: outliers
pvals_outliers <- matrix(nrow=nrow(pooled_chains), ncol=1)
for(i in 1:nrow(pooled_chains)){
  pvals_outliers[i]<- ifelse(simdata[i,1] < obsdata[i,1], 1, 0)
}
mean(pvals_outliers)

# Posterior predictive p-value: skewness
pvals_skew <- matrix(nrow=nrow(pooled_chains), ncol=1)
for(i in 1:nrow(pooled_chains)){
  pvals_skew[i]<- ifelse(simdata[i,2] > obsdata[i,2], 1, 0)
}
mean(pvals_skew)

# Posterior predictive crude pvalue: homoscedasticity
pvals_cor <- matrix(nrow=nrow(pooled_chains), ncol=1)
for(i in 1:nrow(pooled_chains)){
  pvals_cor[i]<- ifelse(simdata[i,3] > obsdata[i,3], 1, 0)
}
mean(pvals_cor)

# Graphical posterior predictive check
plot(obsres, R_obs)
plot(simres, R_obs)

# Trace plots----------------------------------------------------------------------------------------
plot_data <- gather(pooled_chains, key="measure", value="value", c("beta0s", "beta1s", "beta2s", "beta3s", "sigma2s"))

ggplot(plot_data, aes(x = js, y =value)) + 
  geom_line(aes(color = as.factor(chain))) + 
  scale_color_manual(values = c("seagreen4", "orange"),
                     name="Chain")+
  facet_wrap(~measure, scales ="free_y") +
  xlab("") +
  theme_bw()

# Function for autocorrelation----------------------------------------------------------------------------------
acfPlot = function(param, s, plot, title){
  rho<-matrix(0,1) 
  for(i in 0:s){
    rho[i+1] = cor(c(param[1:(1000-i)]), c(param[(i+1):1000]))}
  if(plot){
    df = data.frame(ACF = rho, Lag=0:s)
    graph = ggplot(df, aes(x = Lag, y= ACF)) + 
      geom_bar(width = 0.3, stat = "identity", position = "identity", fill = "seagreen4") + 
      geom_point(color = "orange")+
      labs(title = title) +
      theme_bw()
    plot(graph)}}

p6 <- acfPlot(param = pooled_chains$beta0s,  s = 50, plot = T, title = "beta0")
p7 <- acfPlot(param = pooled_chains$beta1s,  s = 50, plot = T, title = "beta1")
p8 <- acfPlot(param = pooled_chains$beta2s,  s = 50, plot = T, title = "beta2")
p9 <- acfPlot(param = pooled_chains$beta3s,  s = 50, plot = T, title = "beta3")
p10 <- acfPlot(param = pooled_chains$sigma2s, s = 50, plot = T, title = "sigma2")

acfs<-ggarrange(p6,p7,p8,p9,p10, ncol=2, nrow=3)
annotate_figure(acfs, top = text_grob("Autocorrelation", face = "bold", size = 12))


library(kableExtra)
knitr::kable(pooled_sample, caption = "Table 1: Pooled sample statistics for the posterior distributions of the parameters", align = 'l', digits = 6,full_width=F) %>% 
  kable_classic() %>% 
  row_spec(0, bold = TRUE)

# Model 2
mod_2 <- Gibbs_MH_Sampler(y=sleep$duration,x1=sleep$active,x2=sleep$screen,x3=sleep$onset,n=168,mu.00=0,sigma2.00=0.25,mu.10=0,sigma2.10=0.5,mu.30=0,sigma2.30=0.5,A.0=0.00001,B.0=0.00001)

library(bain)
# Hypothesis 1
m1<-lm(duration~active+screen+interaction, data=sleep)
BF_1 <- bain(m1, hypothesis = "active = 0 & screen = 0 & interaction = 0", standardize = T)
print(BF_1, ci = 0.95)
# Hypothesis 2
m2<-lm(duration~active+screen, data=sleep)
BF_2 <- bain(m2, hypothesis = "screen > active", standardize = T)
print(BF_2, ci = 0.95)
# Hypothesis 3
m3<-lm(duration~active+screen+interaction, data=sleep)
BF_3 <- bain(m3, hypothesis = "interaction > screen &interaction > active", standardize = T)
print(BF_3, ci = 0.95)

# Density Plots ---------------------------------------------------------------------------------
p1<-ggplot(pooled_chains, aes(x=beta0s)) +                                          # BETA 0
  geom_histogram(fill = "white", aes(color = as.factor(chain))) + 
  geom_vline(aes(xintercept=mean(beta0s)), color="black", linetype="dashed") +
  scale_color_brewer(palette = "Dark2") +
  ylab("") + xlab("b0") +
  labs(col="Chain") +
  theme_bw()

p2<-ggplot(pooled_chains, aes(x=beta1s)) +                                          # BETA 1
  geom_histogram(fill = "white", aes(color = as.factor(chain))) + 
  geom_vline(aes(xintercept=mean(beta1s)), color="black", linetype="dashed") +
  scale_color_brewer(palette = "Dark2") +
  ylab("") + xlab("b1") +
  labs(col="Chain") +
  theme_bw()

p3<-ggplot(pooled_chains, aes(x=beta2s)) +                                          # BETA 2
  geom_histogram(fill = "white", aes(color = as.factor(chain))) + 
  geom_vline(aes(xintercept=mean(beta2s)), color="black", linetype="dashed") +
  scale_color_brewer(palette = "Dark2") +
  ylab("") + xlab("b2") +
  labs(col="Chain") +
  theme_bw()

p4<-ggplot(pooled_chains, aes(x=beta3s)) +                                          # BETA 3
  geom_histogram(fill = "white", aes(color = as.factor(chain))) + 
  geom_vline(aes(xintercept=mean(beta3s)), color="black", linetype="dashed") +
  scale_color_brewer(palette = "Dark2") +
  labs(col="Chain") +
  ylab("") + xlab("b3") +
  theme_bw()

p5<-ggplot(pooled_chains, aes(x=sigma2s)) +                                          # SIGMA^2
  geom_histogram(fill = "white", aes(color = as.factor(chain))) + 
  geom_vline(aes(xintercept=mean(sigma2s)), color="black", linetype="dashed") +
  scale_color_brewer(palette = "Dark2") +
  labs(col="Chain") +
  ylab("") + xlab("s2") +
  theme_bw()

dens<-ggarrange(p1,p2,p3,p4,p5, ncol=2,nrow=3, legend="bottom", common.legend = T)
annotate_figure(dens, top = text_grob("Marginal posterior densities", face = "bold", size = 12))
