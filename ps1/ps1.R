########## METRICS PSET 1 ##########

#### Question 2 Part E ####
library(tidyverse)
library(MASS)
library(stargazer)
library(kableExtra)

# Function that takes in sigma and rho, simulates data, and outputs results
func <- function(N, sigma, rho) {
    
    # Generate data
    mu    <- c(0,0)
    sigma <- matrix(c(s^2, r*s, r*s, 1), 2)
    errors <- data.frame(mvrnorm(N, mu = mu, Sigma = sigma))
    
    Y1 <- data.frame(errors[,1])
    Y0 <- data.frame(errors[,2])
    
    df <- cbind(Y1, Y0, errors)
    names(df) <- c("Y1","Y0","U1","U0")
    df <- df %>%
        mutate(D = as.integer(ifelse(U1>U0,1,0)),
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
    
    # Save results
    resultnames <- c("\\beta_\\text{OLS}", "ATE", "ATT", "ATU", "E[Y|D=1]-E[Y|D=0]")
    results <- data.frame("Parameter"= resultnames,
                          "Estimate"=c(beta, ATE, ATT, ATU, diff))
    results$Estimate <- round(results$Estimate, digits=3)
    
    return(results)
}


### Step 1: Results for sigma = 2, rho = -0.5, 0, and 0.5

# Setup
N = 10000

# Results
run1 <- func(N, 0.5, -0.5) %>% rename(Estimate_1 = Estimate)
run2 <- func(N, 0.5,  0.0) %>% rename(Estimate_2 = Estimate)
run3 <- func(N, 0.5,  0.5) %>% rename(Estimate_3 = Estimate)

# Combine results
results <- left_join(run1, run2) %>% left_join(.,run3)

# Save results

knitr::kable(results, format="latex",  booktabs = TRUE, escape=FALSE) %>% 
    kable_styling(position = "center",
                  latex_options = "HOLD_position") %>%
    add_header_above(c("($\\rho$,$\\sigma$)" = 1, "(2,-0.5)" = 1, "(2,0)" = 1, "(2,0.5)" = 1),escape=F) %>%
    save_kable(.,"q2d_part1.tex",float = FALSE)
#https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf


# Part 2: vary sigma, fix rho at 0.5



