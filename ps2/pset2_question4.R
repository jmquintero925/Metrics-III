rm(list=ls())
setwd("~/Dropbox/Grad School/Courses/Metrics III - Magne,Jim/Metrics-III/ps2")

library(data.table)
library(ggplot2)
library(stargazer)

## LOAD DATA -----------------------
dt <- fread(file='./IncomeData.csv')

# Exploration
dim(dt) # 1mil rows, 1 column
summary(dt$Y) # median/mean income = 39k and 41k
# max income 142k, min = 2


## PARAMETERS ---------------

delta   <- 2000  # 1/2 size of excluded range
K       <- 20000 # kink location


## PART A: PLOT HISTOGRAM ---------------
ggplot() +
  geom_histogram(data=dt, aes(x=Y), binwidth = delta) +
  scale_y_continuous(limits=c(0,100000), 
                     breaks = scales::breaks_width(50000), 
                     labels = scales::comma) +
  scale_x_continuous(labels=scales::dollar) +
  theme_minimal() + 
  labs(x="Income", y="Count", title="Distribution of Income") 

## PART B: Bunching Probability ---------------
y1 <- K - delta # low bound of excluded range
y2 <- K + delta # upper bound of excluded range

# Non-parametric estimator of f_(y1), f+(y2)
N_f1  <- dt[between(Y,y1 - delta, y1),.N]
N_f2  <- dt[between(Y,y2, y2 + delta),.N]
N_excl <- dt[between(Y,y1,y2),.N]
N_denom <- dt[between(Y,y1 - delta, y2 + delta),.N]

fminus <- N_f1 / (delta * N_denom)
fplus  <- N_f2 / (delta * N_denom)

# Implied P_K
Pk <- (N_excl / N_denom) - delta * (fminus + fplus)

## PART C: Taxable Elasticity ---------------

