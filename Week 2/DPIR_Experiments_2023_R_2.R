############# Week 2 ##############
## Programmed by: Mats Ahrenshop ##
####### Date: 04 May 2020 #########


###########################
# RANDOMIZATION INFERENCE #
###########################


# To illustrate the RI routine we use the "ri" package developed for Gerber and Green (2012).
# The package is currently archived, hence follow these steps to install it
install.packages("https://cran.r-project.org/src/contrib/Archive/ri/ri_0.9.tar.gz", repo=NULL, type="source")


library(tidyverse)
library(data.table)
library(gtools)
library(lubridate)
library(AER)
library(xtable)
library(pBrackets)
library(Hmisc)
library(car)
library(psych)
library(ggpubr)
library(ri)


set.seed(103648)



#---- PRELIMINARIES ----

## Function to remove compute robust SE
robust <- function(model) {
  X <- model.matrix(model)
  n <- dim(X)[1]
  k <- dim(X)[2]
  u <- matrix(resid(model))
  meat1 <- t(X) %*% diag(diag(crossprod(t(u)))) %*% X
  dfc <- n/(n-k) 
  se <- sqrt(dfc*diag(solve(crossprod(X)) %*% meat1 %*% solve(crossprod(X))))
  return(se)
}



#---- RI WITH TOY DATA ----

# Steps:
# - Simulating all possible random assignments
# - Compute sampling distribution of the test statistic under the null hypothesis
# - calculate p-values by comparing the observed test statistic to the
#   distribution of test statistics under the null hypothesis.


## Let's start with a quick example with simulated data

Y <- sample(x = 1:10, size = 18, replace = TRUE)
Z <- sample(x = 0:1, size = 18, replace = TRUE)

# Enumerate all possible permutations given assignment vector (permutation matrix)
perms <- genperms(Z)

# Compute probability of treatment
probs <- genprobexact(Z)

# Estimate observed ATE
ate <- estate(Y, Z, prob = probs)

# Generate full schedule of potential outcomes under sharp H0
Ys <- genouts(Y, Z, ate = 0)

# Generate randomization distribution
distout <- gendist(Ys, perms, prob = probs)

p <- mean(abs(distout) >= abs(ate))
dispdist(distout, ate)


## Simulation alternative:

N <- 50
m <- 25

d <- ifelse(1:N %in% sample(1:N, m), 1, 0)
Y0 <- runif(N, 0, 1)
Y1 <- Y0 + rnorm(N, 2, 2)
Y <- Y1*d + Y0*(1-d)

cbind(Y0, Y1, d, Y)	# look at your data


## Conduct analysis of actual experiment
## Estimate the ATE

# nonparametric
mean(Y[d == 1]) - mean(Y[d == 0])

# or fitting data to ols model
lm(Y ~ d)

Z <- d

# Compute probability of treatment
probs <- genprobexact(Z)

# Estimate observed ATE
ate <- estate(Y, Z, prob = probs)
  
perms <- genperms(Z, maxiter = 10000)

# Create potential outcomes UNDER THE SHARP NULL H
Ys <- genouts(Y, Z, ate = 0)
  
# Generate the sampling distribution based on schedule of potential outcome
# implied by the sharp null hypothesis
distout <- gendist(Ys, perms, prob = probs)
  
sum(distout >= ate)                 # one-tailed comparison used to calculate p-value (greater than)
sum(abs(distout) >= abs(ate))       # two-tailed comparison used to calculate p-value
  
dispdist(distout, ate)               # display p-values, 95% confidence interval, standard error under the null, and graph the sampling distribution under the null


# Estimation of confidence intervals assuming ATE = estimated ATE

# Create potential outcomes UNDER THE ASSUMPTION THAT ATE=ESTIMATED ATE
Ys <- genouts(Y, Z, ate = ate)

# Generate the sampling distribution  based on the schedule of potential outcomes
# implied by the null hypothesis
distout <- gendist(Ys, perms, prob = probs) 

# Display p-values, 95% confidence interval, standard error under the null, 
# and graph the sampling distribution under the null
dispdist(distout, ate)               



#---- YOUNG 2019 ----

## Load data
dat <- read.csv("https://raw.githubusercontent.com/rayduch/DPIR-Experimental-Methods-2023/main/Week%201/young19.csv")

## Create subset of actual wristband days
datw <- subset(dat, dat$wristband_real == 1)

## models include weights as 1/p and heteroskedastic robust SEs (treatment groups are of unequal size)

## ATEs and robust SEs
mod1g <- lm(prob_act_st ~ treat_assign, data = dat[dat$treat_assign != "TP", ], 
            weights = dat$TG_inv[dat$treat_assign != "TP"])
mod1g$se <- robust(mod1g)

mod1p <- lm(prob_act_st ~ treat_assign, data = dat[dat$treat_assign != "TG", ], 
            weights = dat$TP_inv[dat$treat_assign != "TG"])
mod1p$se <- robust(mod1p)

mod2g <- lm(wristband ~ treat_assign, data = datw[datw$treat_assign != "TP", ], 
            weights = datw$TG_inv[datw$treat_assign != "TP"])
mod2g$se <- robust(mod2g)

mod2p <- lm(wristband ~ treat_assign, data = datw[datw$treat_assign != "TG", ], 
            weights = datw$TP_inv[datw$treat_assign != "TG"])
mod2p$se <- robust(mod2p)


## ATE estimates differ slightly from the ones reported in the paper;
## in paper estimated as part of RI routine,
## but SEs are exactly reproduced

## Now let's quickly revisit results in Young 2019 in light of randomization inference

## Hypotheical -- general fear
t <- subset(dat, complete.cases(dat[, 'prob_act_st']) & dat$treat_assign %in% c('C', 'TG'))

# Extract nrow
mod1g$n <- dim(t)[1]

Z <- t$treat_all
Y <- t$prob_act_st
block <- t$block 

# Generate probabilities of treatment by block
probs <- genprobexact(Z = Z, blockvar = block)

# Calculate ate
mod1g$ate <- estate(Y = Y, Z = Z, prob = probs)

# Enumerate all possible ways of random assignment (permutations)
perms <- genperms(Z, maxiter = 100000)

# Generate schedule of potential outcomes under exact H0 (ate = 0)
Ys <- genouts(Y = Y, Z = Z, ate = 0)

# Generate distribution of ATE's (100,000)
distout <- gendist(Ys, perms)

# p-value 2-sided
mod1g$p <- mean(abs(distout) >= abs(mod1g$ate))

dispdist(distout, mod1g$ate)

## Behavioural -- general fear
t <- subset(datw, complete.cases(datw[, 'wristband']) & datw$treat_assign %in% c('C', 'TG'))
mod2g$n <- dim(t)[1]
Z <- t$treat_all
Y <- t$wristband
block <- t$block
probs <- genprobexact(Z, blockvar = block)
mod2g$ate <- estate(Y, Z, prob = probs)
perms <- genperms(Z, maxiter = 100000)
Ys <- genouts(Y, Z, ate = 0)
distout <- gendist(Ys, perms)
mod2g$p <- mean(abs(distout) >= abs(mod2g$ate))




#########################
# Random Cluster Design #
#########################

rm(list = ls())

# Define a treatment effect
treatment_effect <- 1

# Define the individual ids (i)
person <- 1:10

# Define the cluster indicator (j)
hair_color <- c(rep("black",5),rep("brown",5))

# Define the control outcome (Y0)
outcome_if_untreated <- rnorm(n = 10)

# Define the treatment outcome (Y1)
outcome_if_treated <- outcome_if_untreated + treatment_effect

# Version 1 - Not cluster randomized
# Generate all possible non-clustered assignments of treatment (Z)
non_clustered_assignments <- combn(x = unique(person), m = 5)

# Estimate the treatment effect
treatment_effects_V1 <-
  apply(
    X = non_clustered_assignments,
    MARGIN = 2,
    FUN = function(assignment) {
      treated_outcomes <- outcome_if_treated[person %in% assignment]
      untreated_outcomes <- outcome_if_untreated[!person %in% assignment]
      mean(treated_outcomes) - mean(untreated_outcomes)
    })

# Estimate the true standard error
standard_error_V1 <- sd(treatment_effects_V1)

# Plot the histogram of all possible estimates of the treatment effect
hist(treatment_effects_V1,xlim = c(-1,2.5), breaks = 20)


### Cluster
# Version 2 - Cluster randomized
# Generate all possible assignments of treatment when clustering by hair color (Z)
clustered_assignments <- combn(x = unique(hair_color), m = 1)

# Estimate the treatment effect
treatment_effects_V2 <-
  sapply(
    X = clustered_assignments,
    FUN = function(assignment) {
      treated_outcomes   <- outcome_if_treated[person %in% person[hair_color==assignment]]
      untreated_outcomes <- outcome_if_untreated[person %in% person[!hair_color==assignment]]
      mean(treated_outcomes) - mean(untreated_outcomes)
    }
  )

# Estimate the true standard error
standard_error_V2 <- sd(treatment_effects_V2)

# Plot the histogram of all possible estimates of the treatment effect
hist(treatment_effects_V2, xlim = c(-1,2.5), breaks = 20)
