#Cost-effectiveness HPV vaccines with Bayesian synthesis framework
#Gender-neutral vaccination at age 10

set.seed(1)
library(ggplot2)
library(cowplot)
library(reshape2)
library(gtools) #for rdirichlet function

#Need to load these functions once
source('HPVcompare_functions.R')
source('HPVcompare_RiskReduction_functions.R')

#Simulation of parameters
nr.sim <- 1000
vacc.up <- c(0.5,0.5)
source('HPVcompare_ParametersData.R')

#Base-case scenario (Including herd-immunity, and both AW and RRP)
scenario <- c(2,1,3)
source('HPVcompare_RiskReductions.R')
source('HPVcompare_Computation.R')

#ICER 9v vs 2v with 95% credible interval
round(median(Results$ICERs$'9v_vs_2v.mid'))
round(quantile(Results$ICERs$'9v_vs_2v.mid',0.025))
round(quantile(Results$ICERs$'9v_vs_2v.mid',0.975))
