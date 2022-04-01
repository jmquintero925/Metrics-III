##### METRICS PSET 1 #####

## Question 2 Part E ##
library(tidyverse)
library(MASS)

# Setup
s = 2
r = 0.5
N = 10000

mu    <- c(0,0)
sigma <- matrix(c(s^2, r*s, r*s, 1), 2)

errors <- data.frame(mvrnorm(N, mu = mu, Sigma = sigma))

Y1 <- data.frame(errors[,1])
Y0 <- data.frame(errors[,2])

df <- cbind(Y1, Y0, errors)
names(df) <- c("Y1","Y0","U1","U0")
df <- df %>%
    mutate(D = ifelse(U1>U0,1,0),
           Y = D*Y1 + (1-D)*Y0)

# Regression of Y on D
ols <- lm(Y ~ D, df)
beta <- summary(ols)$coefficients[2,1] #1.38999

# Treatment Effects
df <- df %>% mutate(TE=Y1-Y0)
ATE <- mean(df$TE)                     #-0.009993245
ATT <- mean(df[df$D==1,'TE'])          # 1.416414
ATU <- mean(df[df$D==0,'TE'])          #-1.388713

# Compute parameter
diff <- mean(df[df$D==1,'Y']) - mean(df[df$D==0,'Y'])
sb1  <- diff - ATT
sb2  <- diff - ATU
####Why are they diff selection biases?? Well because we're decomposing it into diff things!

