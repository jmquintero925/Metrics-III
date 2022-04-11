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
RHO1    <- 1.0   # tax rate before kink
RHO2    <- 0.8   # tax rate after kink


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
# N_denom <- dt[between(Y,y1 - delta, y2 + delta),.N]
N_denom <- dt[,.N]

FMINUS <- N_f1 / (delta * N_denom)
FPLUS  <- N_f2 / (delta * N_denom)

# Implied P_K
Pk <- (N_excl / N_denom) - delta * (FMINUS + FPLUS)

## PART C: Taxable Elasticity ---------------
PARAM <- c(FMINUS,FPLUS,RHO1,RHO2,K, y1, y2)

# Three specifications
compute_integral <- function(beta, m, param=PARAM){
  
  # unpack parameters
  fminus <- param[1]
  fplus <- param[2]
  rho1 <- param[3]
  rho2 <- param[4]
  K <- param[5]
  y1 <- param[6]
  y2 <- param[7]
  
  # compute eta1, eta2
  eta1 <- y1 * (rho1 ^ (-beta))
  eta2 <- y2 * (rho2 ^ (-beta))
  
  # compute boundary densities
  phi1 <- fminus * (rho1^beta)
  phi2 <- fplus  * (rho2^beta)
  
  # compute integral
  if ( m==1) {
    ans = (eta2-eta1) * ((1/2)*phi1 + (1/2)*phi2)
  } else if ( m==2) {
    ans = (eta2-eta1) * ((7/8)*phi1 + (1/8)*phi2)
  } else if ( m==3) {
    ans = (eta2-eta1) * ((1/8)*phi1 + (7/8)*phi2)
  }
  return(ans)
}

# PLOT
betagrid <- seq(from=0,to=2,length.out=100)
m1 = compute_integral(betagrid, m=1)
m2 = compute_integral(betagrid, m=2)
m3 = compute_integral(betagrid, m=3)

dt_plot <- data.table(
  beta = betagrid,
  m1 = m1,
  m2 = m2,
  m3 = m3
)

dt_plot <- melt(dt_plot, id.vars='beta') # reshape long


ggplot(dt_plot) +
  geom_point(aes(x=beta,y=value, color=variable)) +
  # geom_point(aes(x=beta,y=m2), color='blue') +
  # geom_point(aes(x=beta,y=m3), color='black') +
  geom_hline(yintercept = Pk) +
  scale_color_manual(values = c("red", "blue", 'black'), name="Specification")+
  theme_minimal(base_size = 15) +
  labs(x="Beta", y="PK", title="Relationship between Kink and Elasticity") +
  theme(plot.title.position = "plot", 
        legend.position = "bottom")









