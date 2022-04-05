########## METRICS PSET 1 ##########
setwd("/Users/Isaac/Documents/Code/Metrics III/ps1/")

#### Question 2 Part E ####
library(tidyverse)
library(MASS)
library(stargazer)
library(kableExtra)
library(latex2exp)
library(sjmisc)
set.seed(1126)

# Function that takes in sigma and rho, simulates data, and outputs results
func <- function(N, s, r) {
    
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
    beta <- summary(ols)$coefficients[2,1]
    
    # Treatment Effects
    df <- df %>% mutate(TE=Y1-Y0)
    ATE <- mean(df$TE)
    ATT <- mean(df[df$D==1,'TE'])
    ATU <- mean(df[df$D==0,'TE'])
    
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
s <- 2

# Results
run1 <- func(N, s, -0.5) %>% rename(Estimate_1 = Estimate)
run2 <- func(N, s,  0.0) %>% rename(Estimate_2 = Estimate)
run3 <- func(N, s,  0.5) %>% rename(Estimate_3 = Estimate)

# Combine results
results <- left_join(run1, run2) %>% left_join(.,run3)
names(results) <- c("($\\sigma$,$\\rho$)", "(2,-0.5)", "(2,0)", "(2,0.5)")

# Save results
knitr::kable(results, format="latex",  booktabs = TRUE, escape=FALSE) %>% 
    kable_styling(position = "center",
                  latex_options = "HOLD_position") %>%
    add_header_above(c("Parameter"=1,"Estimate"=3),escape=F) %>%
    save_kable(.,"./Tables/q2d_tab1.tex",float = FALSE)

# Part 2: vary sigma, fix rho at 0.5
rho = 0.5
sigmas <- seq(from = 0, to = 5, by = 0.01)

resultnames <- c("\\beta_\\text{OLS}", "ATE", "ATT", "ATU", "E[Y|D=1]-E[Y|D=0]")
df2 <- data.frame("Parameter"= resultnames)
for (i in 1:length(sigmas)) {
    newcolname <- paste("Est_",as.character(i),sep="")
    run <- func(N, sigmas[i], -0.5) %>% rename(!!newcolname := Estimate)
    df2 <- suppressMessages(left_join(df2, run))
}

df2 <- sjmisc::rotate_df(head(df2),, cn=TRUE) %>% 
    rename(beta = `\\beta_\\text{OLS}`,
           param = `E[Y|D=1]-E[Y|D=0]`)
df2$sigma <- sigmas

betaplot <- ggplot(data=df2, aes(x=sigma, y=beta)) + 
    geom_point() + 
    ylab(TeX(r'(ATU)')) +
    xlab(TeX(r'($\sigma$)'))
ggsave(betaplot, filename="./Figures/p2_qe_beta.pdf")

ATEplot <-  ggplot(data=df2, aes(x=sigma, y=ATE)) + 
    geom_point() + 
    ylab(TeX(r'(ATU)')) +
    xlab(TeX(r'($\sigma$)'))
ggsave(ATEplot, filename="./Figures/p2_qe_ATE.pdf")

ATTplot<- ggplot(data=df2, aes(x=sigma, y=ATT)) + 
    geom_point() + 
    ylab(TeX(r'(ATU)')) +
    xlab(TeX(r'($\sigma$)'))
ggsave(ATTplot, filename="./Figures/p2_qe_ATT.pdf")

ATUplot<- ggplot(data=df2, aes(x=sigma, y=ATU)) + 
    geom_point() + 
    ylab(TeX(r'(ATU)')) +
    xlab(TeX(r'($\sigma$)'))
ggsave(ATUplot, filename="./Figures/p2_qe_ATU.pdf")


