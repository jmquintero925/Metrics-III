rm(list=ls())
setwd("~/Dropbox/Grad School/Research/3_TestOptional/2_Code/src")

library(data.table)
library(ggplot2)
library(stargazer)

## LOAD DATA -----------------------
load(file='./clean_scorecard.Rdata')
load(file='../output/test_optional_policy.Rdata')


## Merge Scorecard and Test Optional DATA -----------------------
dt <- merge(dt,TOdata, by='UNITID', all.x=TRUE)
  

## CLEAN AND PREP -----------------------

# Define treatment vars
dt[, year_since_treatment := year - year_switch] # relative time

# Categorize
dt[, treated := if_else(year_switch < 2020, 'treated pre covid','treated 2020')]
dt[is.na(year_switch), 'treated'] = 'never treated'

## DESCRIPTIVE STATS -----------------------
aux_mean <- function(x) mean(x, na.rm=TRUE) # drop na mean function

# COMPARE TREATED VS. NEVER TREATED VS. TREATED 2020
dt_summarystats <- dt %>%
  filter(year==2011) %>%
  group_by(treated) %>%
  mutate(N=n()) %>% 
  summarize_all(aux_mean)

dt_summarystats <- data.frame(t(dt_summarystats))
dt_summarystats <- drop_na(dt_summarystats)

dt_summarystats <- dt_summarystats %>%
  rename(`Test Required` = X1, `Test optional during 2020` = X2, `Test Optional prior to 2020` = X3) 

stargazer(dt_summarystats, type='latex', digits = 2,
          title='Summary statistics for 2011',
          out='../output/tabs/descriptive.tex', 
          summary=FALSE)


# COMPARE COHORTS WITHIN TREATED < 2020
dt_summarystats <- dt[treated=='treated pre covid'] # only schools TO pre covid

dt_summarystats <- dt_summarystats %>%
  filter(year==2011) %>%
  filter(between(year_switch, 2011, 2019)) %>%
  group_by(year_switch) %>%
  mutate(N=n()) %>% 
  summarize_all(aux_mean) %>%
  select_if(~ !any(is.na(.)))

dt_summarystats <- data.frame(t(dt_summarystats)) # transpose
dt_summarystats <- drop_na(dt_summarystats)

stargazer(dt_summarystats, type='latex', digits = 2,
          title='Summary statistics for 2011 by Cohort',
          out='../output/tabs/descriptive_cohort.tex', summary=FALSE)


## SAMPLE SELECTION -----------------------
save(dt, file='../output/full.Rdata')
dt <- dt[treated=='treated pre covid'] # only schools TO pre covid
dt <- dt[between(year,2011,2019)] # restrict 2011-2019
dt <- dt[!is.na(adm_rate)]
dt <- dt[!is.na(avg_cost)]

## Save out cleaned data -----------------------
save(dt, file='../output/main.Rdata')
