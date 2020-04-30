library(foreign)
library(lme4)
library(pbkrtest)

d<-as.data.frame(read.spss('study1.sav'))
m0 <- lmer(m ~ Condition + (1 | PPNR), data = d)
summary(m0)

#EXP 1
condition <- c(rep(-.5, 400), rep(.5, 400)) 
PPNR <- rep(c(1:20), 40)
X <- as.data.frame(cbind(PPNR = PPNR, condition = condition))
X
b <- c(2, .1) # fixed intercept and slope
V1 <- .9 # random intercept variance
s <- .11 # residual variance
model1 <- makeLmer(dv ~ condition + (1|PPNR), fixef=b, VarCorr=V1, sigma=s, data=X)
summary(model1)


model1 <- makeLmer(y ~ x + (1|g), fixef=b, VarCorr=V1, sigma=s, data=X)
model2 <- makeGlmer(z ~ x + (x|g), family="poisson", fixef=b, VarCorr=V2, data=X)
powerSim(model1, nsim=20)
powerSim(model2, nsim=20)

##########
Subject <- gl(n = 20, k = 10)
C1<-c(-.5,.5)
C2<-c(-.5,.5)
# creates "frame" for our data
X <- expand.grid(Subject=Subject, C1=C1, C2=C2)

b3 <- c(27, 1, 1,.4) # fixed intercept and slope
SubVC3  <- matrix(c(17.2,.73, .73, 0), 2) # random intercept and slope variance-covariance matrix
s3 <- .25 # residual sd 

mod <- makeLmer(DV ~ C1 * C2 + (1|Subject), 
                fixef=b3, VarCorr = SubVC3, sigma=s3, data=X)
summary(mod)


mod0 <- makeLmer(DV ~ C1 + C2 + (1|Subject), 
                 fixef=b3[-4], VarCorr = SubVC3, sigma=s3, data=X)


simPwr <- powerSim(mod, test = compare(mod0, method="lr"), nsim = 100, seed = 1)
simPwr

######################

library("lme4")        # model specification / estimation
library("afex")        # anova and deriving p-values from lmer
library("broom.mixed") # extracting data from model fits 
library("faux")        # data simulation
# NOTE: to install the 'faux' package, use:
devtools::install_github("debruine/faux")
library("tidyverse")   # data wrangling and visualisation

set.seed(123)
# set up the custom data simulation function
my_sim_data <- function(
  nsubj  = 200, # number of subjects
  nitem  = c(inclusion = 24, exclusion = 24),  # number of items
  b0     = 27.50, # the average internal temperature is 37.0 °C
  b1     =  -0.3, # effect of Exclusion on temperature
  I0i_sd =  0, # by-item random intercept sd
  S0s_sd = 0.3, # by-subject random intercept sd  random variation in 
  #temperature for each subject
  S1s_sd =  0, # by-subject random slope sd
  scor   = 0, # correlation between intercept and slope
  err_sd = 0.6){  # residual (standard deviation)  twice the size 
  #of the by-subject random intercept (S0s_sd) SD by default 
  
  # simulate items
  items <- faux::sim_design(
    between = list(category = c("inclusion", "exclusion")),
    n = nitem,
    sd = I0i_sd,
    dv = "I0i",
    id = "time",
    plot = FALSE
  )
  
  # effect code category
  items$cat <- recode(items$category, "inclusion" = 0.5, "exclusion" = -0.5)
  
  # simulate subjects
  subjects <- faux::sim_design(
    within = list(effect = c(S0s = "By-subject random intercepts", 
                             S1s = "By-subject random slopes")), 
    n = nsubj,
    sd = c(S0s_sd, S1s_sd), 
    r = scor,
    id = "subj_id",
    plot = FALSE
  )
  
  # simulate trials

  dat_sim <- crossing(subj_id = subjects$subj_id,
                      time = items$time) %>%
    inner_join(subjects, "subj_id") %>%
    inner_join(items, "time") %>%
    mutate(err = rnorm(nrow(.), mean = 0, sd = err_sd)) %>%
    mutate(RT = b0 + I0i + S0s + (b1 + S1s) * cat + err) %>%
    select(subj_id, time, category, cat, RT)
  
  dat_sim
}

data <- my_sim_data() 
View(data)

# Load libraries
#+message = FALSE
if (!require(simr)) {install.packages('simr')}
if (!require(lme4)) {install.packages('lme4')}
if (!require(lmerTest)) {devtools::install_github('cran/lmerTest')}
if (!require(r2glmm)) {install.packages('r2glmm')}
if (!require(effects)) {install.packages('effects')}

model0 <- lmer(RT ~ + (1|subj_id), data = data)
model <- lmer(RT ~ cat + (1|subj_id), data = data)
summary(model)

# Set the smallest effect size of interest (SESOI)
SESOI <- -.025 #SESOI 
fixef(model)["cat"] <- SESOI

# POWER - Main effect
# Power to detect SESOI 
pwr1 <- powerSim(model, test = compare(model0, method = "lr"), nsim = 15, progress = T)
nsim
pwr1

powerSim
# Power curve, keeping SESOI, for varying N
pwr2 <- extend(model, along = "subj_id", n = 150)
power.curve <- powerCurve(pwr2, test = compare(model0, method = "lr"), nsim = 15, along = "subj_id", progress = T)
power.curve
plot(power.curve)
