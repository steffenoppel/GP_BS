#### BREEDING SUCCESS OF GIANT PETRELS



# setting working directory
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/Students/MichelleRisi/GP_BS"), silent=T)

# 1. LOADING PACKAGES ------
library(tidyverse)
library(dplyr)
library(lubridate)
library(data.table)
library(dtplyr)
library(climwin)
library(betareg)
select<-dplyr::select
filte<-dplyr::filter


# 2. LOADING DATA ---------
GP<-fread("data/Master_Breeding_Success_2000_2023_all islands.csv") %>%
  gather(key="SPECLOC", value="BreedSuccess",-Year) %>%
  separate_wider_delim(cols=SPECLOC,delim="_", names=c("Species","Island")) %>%
  filter(!is.na(BreedSuccess))
  


# 3. SIMPLE PLOT OF DATA ---------

ggplot(data=GP,aes(x=Year,y=BreedSuccess, col=Island)) +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  facet_wrap(~Species, ncol=1)


# 4. SIMPLE CORRELATION MATRIX ---------

GPmat<-as.matrix(fread("data/Master_Breeding_Success_2000_2023_all islands.csv") %>%
  select(-Year))

cor(GPmat, use = "pairwise.complete.obs") ## overall fairly weak correlations, SGPE Gough-Marion is the highest



# 5. FIT SIMPLE LINEAR MIXED MODEL WITH ISLAND; SPECIES AND YEAR EFFECTS ---------
# beta regression most appropriate for proportions
# https://topmodels.r-forge.r-project.org/betareg/vignettes/betareg.html


m0<-betareg(BreedSuccess ~Species + Island + as.factor(Year), data=GP, link="logit")
summary(m0)

## you need to add a sensible environmental covariate
GP$SST<-runif(dim(GP)[1],3,7)
m1<-betareg(BreedSuccess ~Species + Island + as.factor(Year) + SST, data=GP, link="logit")
summary(m1)

## compare models
AIC(m1)
AIC(m0) ## currently better (lower AIC) because the SST variable is randomly generated and explains nothing





# 6. DOWNLOAD OCEAN DATA ---------
# first download data: https://theoceancode.netlify.app/post/dl_env_data_r/
# or here: https://cran.r-project.org/web/packages/oceanexplorer/vignettes/oceanexplorer.html




# 7. EXTRACT A SENSIBLE WINDOW OF CLIMATE DATA ---------
# first fit a basic model as above
# then find best time window to explain breeding success: 

wing_50_basemod <- lme4::lmer(
  wing_50 ~
    brood_size + hatch_doy_sc +
    (1 | year_f) + (1 | nestcode_rearing),
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+06)),
  data = data_nestlings_wing_50
)


wing_50_sw1 <- slidingwin(
  baseline = wing_50_basemod,
  xvar = list(
    temp_mean = data_weather$t_daily_mean
  ),
  stat = c("mean"),
  func = c("lin", "quad"),
  type = "relative", # the window is relative to the individual
  range = c(50, 0),
  # 0 is the date given in bdate, 50 is the number of days before that date
  cinterval = "day",
  cdate = data_weather$date,
  bdate = data_nestlings_wing_50$date_50days
)




