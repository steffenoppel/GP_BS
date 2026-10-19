###########################################################################
##### GIANT PETREL WINTER ATTENDANCE ANALYSIS --------  ############################
###########################################################################

## attempt to create model for Northern Giant Petrels on Marion Island relating winter attendance to breeding performance
## work started by Michelle Risi
## model conceived by Steffen Oppel
## two components - multi-state model for breeding season transitions and JS model for winter attendance


# Clear workspace ---------------------------------------------------------

rm(list=ls()) 
Sys.setenv(LANG = "en") ##change language of error messages to english

# Load libraries, functions etc -------------------------------------------
library(lubridate)
library(tidyverse)
library(dplyr)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select
library(coda)
library(MCMCvis)
library(tictoc)
library(nimble)
library(basicMCMCplots) # for trace plots called chainsPlot


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/General - VOWA_PROJ_Golden_Eagle/ANALYSES/GOEA_tradeoff"), silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/General - VOWA_PROJ_Golden_Eagle/ANALYSES/GOEA_tradeoff"), silent=T)


# Load prepared data ---------------------------------------------------------------





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP NIMBLE CODE FOR ANNUAL BREEDING TRANSITION MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specify model in NIMBLE format
GP.breed.model<-nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi: yearly survival probability
  # p.breed: probability to breed
  # p.succ: probability to raise a chick given that breeding occurred
  
  # p.obs: probability to be observed at all
  # p.obs.breed: probability to be observed as breeder
  # p.obs.succ: probability to be observed as successful breeder
  
  # -------------------------------------------------
  # States (S):
  # 1 not breeding
  # 2 failed breeding
  # 3 successful breeding
  # 4 dead
  
  # Observations (O):
  # 1 observed as loafer
  # 2 observed as failed breeder
  # 3 observed as successful breeder
  # 4 not observed
  
  # -------------------------------------------------
  
  # Priors and constraints
  
  #### SLOPE PARAMETERS FOR LINEAR PREDICTORS
  b.phi.1 ~ dnorm(0, tau=0.01)           # Prior for dependence effect on survival probability on logit scale
  b.breed.win ~ dnorm(0, tau=0.01)        # Prior for season effect on survival probability on logit scale
  b.succ.win ~ dnorm(0, tau=0.01)         # Prior for PFDP effect on survival probability on logit scale
  b.obs.eff ~ dnorm(0, tau=0.01)         # Prior for effort effect on observing bird at all
  b.obs.breed.eff ~ dnorm(0, tau=0.01)         # Prior for effort effect on observing bird breeding
  b.obs.succ.eff ~ dnorm(-2, tau=0.01)         # Prior for effort effect on observing bird successfully breeding
  
  #### INTERCEPTS FOR LINEAR PREDICTORS
  for(ag in 1:3) {
    mean.phi[ag] ~ dbeta(15, 2)   # weakly informed prior intercept for annual survival probabilities
    lp.mean[ag] <- log(mean.phi[ag]/(1 - mean.phi[ag]))    # logit transformed survival intercept
  }
  
  mean.breed ~ dbeta(2, 2)   # fairly uninformative prior intercept for annual breeding probabilities
  lp.breed <- log(mean.breed/(1 - mean.breed))    # logit transformed survival intercept
  
  mean.succ ~ dbeta(2, 2)   # fairly uninformative prior intercept for annual success probabilities
  lp.succ <- log(mean.succ/(1 - mean.succ))    # logit transformed survival intercept
  
  p.obs.int ~ dbeta(2,2)     # uninformative prior intercept for all observation probabilities
  lp.p.obs.int <- log(p.obs.int/(1 - p.obs.int))    # logit transformed observation probability intercept
  
  p.obs.breed.int ~ dbeta(2,2)     # uninformative prior intercept for all breeding observation probabilities
  lp.p.obs.breed <- log(p.obs.breed.int/(1 - p.obs.breed.int))    # logit transformed observation probability intercept
  
  p.obs.succ.int ~ dbeta(2,2)     # uninformative prior intercept for all success observation probabilities
  lp.p.obs.succ <- log(p.obs.succ.int/(1 - p.obs.succ.int))    # logit transformed observation probability intercept
  
  
  #### RANDOM EFFECTS 
  for (t in 1:(n.years)){
    ran.eff.year[t] ~ dnorm(0, sd=0.2)						# random variation around annual survival
  }
  

  #### STATE TRANSITION PROBABILITIES
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      
      ## survival 
      logit.phi[i,t] <- lp.mean[age[i,t]] +       ### intercept for mean survival (vary by age, sex, or breeding experience)
        ran.eff.year[year[t]]                     ### random year effect
      phi[i,t]<-ilogit(logit.phi[i,t])
      
      ## breeding 
      logit.p.breed[i,t] <- lp.breed + b.breed.win*win.att[i,t]      ### linear predictor for prob of breeding including a winter attendance effect
      p.breed[i,t]<-ilogit(logit.p.breed[i,t])
      
      ## success 
      logit.p.succ[i,t] <- lp.succ + b.succ.win*win.att[i,t]      ### linear predictor for prob of success including a winter attendance effect
      p.succ[i,t]<-ilogit(logit.p.succ[i,t])
      
    } #t
  } #i
  


  #### OBSERVATION PROBABILITY

  for (i in 1:nind){
    for (t in f[i]:(n.occasions)){

      logit.p.obs[i,t] <- lp.p.obs.int + b.obs.eff*EFFORT[i,t] ## NEED TO COMPLETE A SENSIBLE OBSERVATION MODEL
      p.obs[i,t]<-ilogit(logit.p.obs[i,t])
      
      logit.p.obs.breed[i,t] <-lp.p.obs.breed + b.obs.breed.eff*EFFORT[i,t]  ## NEED TO COMPLETE A SENSIBLE BREEDING OBSERVATION MODEL
      p.obs.breed[i,t]<-ilogit(logit.p.obs.breed[i,t])
      
      logit.p.obs.succ[i,t] <- lp.p.obs.succ + b.obs.succ.eff*EFFORT[i,t]   ## NEED TO COMPLETE A SENSIBLE SUCCESS OBSERVATION MODEL
      p.obs.succ[i,t]<-ilogit(logit.p.obs.succ[i,t])
      
    } #t
  } #i
  
  
  # -------------------------------------------------
  # Define state-transition and observation matrices 
  # -------------------------------------------------
  
  for (i in 1:nind){
    
    for (t in f[i]:(n.occasions-1)){
      
      # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
      
      ps[1,i,t,1]<-phi[i,t] * (1-p.breed[i,t])  
      ps[1,i,t,2]<-phi[i,t] * p.breed[i,t] * (1-p.succ[i,t]) 
      ps[1,i,t,3]<-phi[i,t] * p.breed[i,t] * p.succ[i,t]
      ps[1,i,t,4]<-(1-phi[i,t])
      
      ps[2,i,t,1]<-phi[i,t] * (1-p.breed[i,t]) 
      ps[2,i,t,2]<-phi[i,t] * p.breed[i,t] * (1-p.succ[i,t]) 
      ps[2,i,t,3]<-phi[i,t] * p.breed[i,t] * p.succ[i,t]
      ps[2,i,t,4]<-(1-phi[i,t])
      
      ps[3,i,t,1]<-phi[i,t] * (1-p.breed[i,t])
      ps[3,i,t,2]<-phi[i,t] * p.breed[i,t] * (1-p.succ[i,t])
      ps[3,i,t,3]<-phi[i,t] * p.breed[i,t] * p.succ[i,t]
      ps[3,i,t,4]<-(1-phi[i,t])
      
      ps[4,i,t,1]<-0
      ps[4,i,t,2]<-0
      ps[4,i,t,3]<-0
      ps[4,i,t,4]<-1 ## dead birds stay dead
      
      # Define probabilities of O(t) [last dim] given S(t)  [first dim]
      
      po[1,i,t,1]<-(1-p.obs.breed[i,t])*p.obs[i,t]
      po[1,i,t,2]<-0
      po[1,i,t,3]<-0
      po[1,i,t,4]<-(1-p.obs[i,t])
      
      po[2,i,t,1]<-p.obs[i,t]*(1-p.obs.breed[i,t])
      po[2,i,t,2]<-p.obs[i,t]*p.obs.breed[i,t]*(1-p.obs.succ[i,t])
      po[2,i,t,3]<-0
      po[2,i,t,4]<-(1-p.obs[i,t])
      
      po[3,i,t,1]<-p.obs[i,t]*(1-p.obs.breed[i,t])
      po[3,i,t,2]<-p.obs[i,t]*p.obs.breed[i,t]*(1-p.obs.succ[i,t])
      po[3,i,t,3]<-p.obs[i,t]*p.obs.breed[i,t]*p.obs.succ[i,t]
      po[3,i,t,4]<-(1-p.obs[i,t])
      
      po[4,i,t,1]<-0
      po[4,i,t,2]<-0
      po[4,i,t,3]<-0
      po[4,i,t,4]<-1 ## we omit the probability to have a carcass recovered - just to keep the model simple. Assumption is that dead birds are never seen
      
    } #t
  } #i
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 2 ## alive when first marked
    for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:4])     ### check different CJS distribution in nimbleEcology dcjs
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,1:4])
    } #t
  } #i
  
}) ## end of nimble code chunk













