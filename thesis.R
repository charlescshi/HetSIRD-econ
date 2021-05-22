library(ggplot2)
library(dplyr)
library(numDeriv)
library(dfoptim)
library(lubridate)
library(foreach)
library(purrr)
library(optimParallel)

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))
shell(paste("PowerShell -Command \"& {(Get-Process -id ",Sys.getpid(),").ProcessorAffinity = ",affinity.mask,"}\"",sep=""))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)

## Parameters

hs = 1
hia = 1
his = 1
hih = 0.8
beta = 0.4
gamma_a = 1/6
gamma_s = 1/10
gamma_h = 1/23
sigma = 1
rho = 0.05/365
delta = 0.67/365
tau_a_y = 0.79
tau_s_y = 0.20
tau_h_y = 0.01
tau_a_m = 0.5
tau_s_m = 0.476
tau_h_m = 0.024
tau_a_o = 0.31
tau_s_o = 0.614
tau_h_o = 0.076
pi_a_y = 0
pi_s_y = 0
pi_h_y = 0.02
pi_a_m = 0
pi_s_m = 0
pi_h_m = 0.079
pi_a_o = 0
pi_s_o = 0.123
pi_h_o = 0.188
v = 31755
psi = 1021
kappa_a_y = pi_a_y*v
kappa_s_y = pi_s_y*v
kappa_h_y = pi_h_y*v
kappa_a_m = pi_a_m*v
kappa_s_m = pi_s_m*v
kappa_h_m = pi_h_m*v
kappa_a_o = pi_a_o*v
kappa_s_o = pi_s_o*v
kappa_h_o = pi_h_o*v

## Global settings

horizon = 600
chunk <- function(x, n) (mapply(function(a, b) (x[a:b]),
                                seq.int(from=1, to=length(x), by=n),
                                pmin(seq.int(from=1, to=length(x), by=n)+(n-1), length(x)),
                                SIMPLIFY=FALSE))

## -----------------
## Auxillary functions
## -----------------

u_sus = function(a){
  log(a) - a + hs
}

u_asym = function(a){
  log(a) - a + hia
}

u_sym = function(a){
  log(a) - a + his
}

u_hosp = function(a){
  log(a) - a + hih
}

contact = function(c1, c2){
  c1*c2
}

## -------------------------------------------------------
##
## Optimal path
##
## -------------------------------------------------------

## Disease dynamics

## sird_vector: young (s,ia,is,ih,r,d), mid (s,i..,r,d), old (s,i..,r,d)
## contacts_vector: activity young (s,ia,is,ih,r,d), mid (s,i..,r,d), old (s,i..,r,d)
## this version should be used to chart out the trajectory.
next_period_sird = function(sird, c){
  next_sird = vector(mode="numeric", length=18)
  temp_young = sigma*(contact(c[1], c[2])*sird[2]
                      + contact(c[1], c[3])*sird[3]
                      + contact(c[1], c[4])*sird[4]
                      + contact(c[1], c[8])*sird[8]
                      + contact(c[1], c[9])*sird[9]
                      + contact(c[1], c[10])*sird[10]
                      + contact(c[1], c[14])*sird[14]
                      + contact(c[1], c[15])*sird[15]
                      + contact(c[1], c[16])*sird[16])
  temp_mid = sigma*(contact(c[7], c[2])*sird[2]
                      + contact(c[7], c[3])*sird[3]
                      + contact(c[7], c[4])*sird[4]
                      + contact(c[7], c[8])*sird[8]
                      + contact(c[7], c[9])*sird[9]
                      + contact(c[7], c[10])*sird[10]
                      + contact(c[7], c[14])*sird[14]
                      + contact(c[7], c[15])*sird[15]
                      + contact(c[7], c[16])*sird[16])
  temp_old = sigma*(contact(c[13], c[2])*sird[2]
                      + contact(c[13], c[3])*sird[3]
                      + contact(c[13], c[4])*sird[4]
                      + contact(c[13], c[8])*sird[8]
                      + contact(c[13], c[9])*sird[9]
                      + contact(c[13], c[10])*sird[10]
                      + contact(c[13], c[14])*sird[14]
                      + contact(c[13], c[15])*sird[15]
                      + contact(c[13], c[16])*sird[16])
  next_sird[1] = sird[1]-beta*sird[1]*temp_young
  next_sird[2] = sird[2]+beta*tau_a_y*sird[1]*temp_young - gamma_a*sird[2]
  next_sird[3] = sird[3]+beta*tau_s_y*sird[1]*temp_young - gamma_s*sird[3]
  next_sird[4] = sird[4]+beta*tau_h_y*sird[1]*temp_young - gamma_h*sird[4]
  next_sird[5] = sird[5]+(1-pi_a_y)*gamma_a*sird[2] + (1-pi_s_y)*gamma_s*sird[3] + (1-pi_h_y)*gamma_h*sird[4]
  next_sird[6] = sird[6]+(pi_a_y)*gamma_a*sird[2] + (pi_s_y)*gamma_s*sird[3] + (pi_h_y)*gamma_h*sird[4]
  next_sird[7] = sird[7]-beta*sird[7]*temp_mid
  next_sird[8] = sird[8]+beta*tau_a_m*sird[7]*temp_mid - gamma_a*sird[8]
  next_sird[9] = sird[9]+beta*tau_s_m*sird[7]*temp_mid - gamma_s*sird[9]
  next_sird[10] = sird[10]+beta*tau_h_m*sird[7]*temp_mid - gamma_h*sird[10]
  next_sird[11] = sird[11]+(1-pi_a_m)*gamma_a*sird[8] + (1-pi_s_m)*gamma_s*sird[9] + (1-pi_h_m)*gamma_h*sird[10]
  next_sird[12] = sird[12]+(pi_a_m)*gamma_a*sird[8] + (pi_s_m)*gamma_s*sird[9] + (pi_h_m)*gamma_h*sird[10]
  next_sird[13] = sird[13]-beta*sird[13]*temp_old
  next_sird[14] = sird[14]+beta*tau_a_o*sird[13]*temp_old - gamma_a*sird[14]
  next_sird[15] = sird[15]+beta*tau_s_o*sird[13]*temp_old - gamma_s*sird[15]
  next_sird[16] = sird[16]+beta*tau_h_o*sird[13]*temp_old - gamma_h*sird[16]
  next_sird[17] = sird[17]+(1-pi_a_o)*gamma_a*sird[14] + (1-pi_s_o)*gamma_s*sird[15] + (1-pi_h_o)*gamma_h*sird[16]
  next_sird[18] = sird[18]+(pi_a_o)*gamma_a*sird[14] + (pi_s_o)*gamma_s*sird[15] + (pi_h_o)*gamma_h*sird[16]
  next_sird
}

## sird_vector: young (s,ia,is,ih), mid (s,i..), old (s,i..)
## contacts_vector: activity young (s/ia,is,ih), mid (s,i..), old (s,i..)
## this version should be used during optimiaztion.
next_period_sird_cond = function(sird, c){
  next_sird = vector(mode="numeric", length=12)
  temp_young = sigma*(contact(c[1], c[1])*sird[2]
                      + contact(c[1], c[2])*sird[3]
                      + contact(c[1], c[3])*sird[4]
                      + contact(c[1], c[4])*sird[6]
                      + contact(c[1], c[5])*sird[7]
                      + contact(c[1], c[6])*sird[8]
                      + contact(c[1], c[7])*sird[10]
                      + contact(c[1], c[8])*sird[11]
                      + contact(c[1], c[9])*sird[12])
  temp_mid = sigma*(contact(c[4], c[1])*sird[2]
                    + contact(c[4], c[2])*sird[3]
                    + contact(c[4], c[3])*sird[4]
                    + contact(c[4], c[4])*sird[6]
                    + contact(c[4], c[5])*sird[7]
                    + contact(c[4], c[6])*sird[8]
                    + contact(c[4], c[7])*sird[10]
                    + contact(c[4], c[8])*sird[11]
                    + contact(c[4], c[9])*sird[12])
  temp_old = sigma*(contact(c[7], c[1])*sird[2]
                    + contact(c[7], c[2])*sird[3]
                    + contact(c[7], c[3])*sird[4]
                    + contact(c[7], c[4])*sird[6]
                    + contact(c[7], c[5])*sird[7]
                    + contact(c[7], c[6])*sird[8]
                    + contact(c[7], c[7])*sird[10]
                    + contact(c[7], c[8])*sird[11]
                    + contact(c[7], c[9])*sird[12])
  next_sird[1] = sird[1]-beta*sird[1]*temp_young
  next_sird[2] = sird[2]+beta*tau_a_y*sird[1]*temp_young - gamma_a*sird[2]
  next_sird[3] = sird[3]+beta*tau_s_y*sird[1]*temp_young - gamma_s*sird[3]
  next_sird[4] = sird[4]+beta*tau_h_y*sird[1]*temp_young - gamma_h*sird[4]
  next_sird[5] = sird[5]-beta*sird[5]*temp_mid
  next_sird[6] = sird[6]+beta*tau_a_m*sird[5]*temp_mid - gamma_a*sird[6]
  next_sird[7] = sird[7]+beta*tau_s_m*sird[5]*temp_mid - gamma_s*sird[7]
  next_sird[8] = sird[8]+beta*tau_h_m*sird[5]*temp_mid - gamma_h*sird[8]
  next_sird[9] = sird[9]-beta*sird[9]*temp_old
  next_sird[10] = sird[10]+beta*tau_a_o*sird[9]*temp_old - gamma_a*sird[10]
  next_sird[11] = sird[11]+beta*tau_s_o*sird[9]*temp_old - gamma_s*sird[11]
  next_sird[12] = sird[12]+beta*tau_h_o*sird[9]*temp_old - gamma_h*sird[12]
  next_sird
}

## Objective in each time period
## sird_vector: young (s,ia,is,ih), mid (s,i..), old (s,i..)
## contacts_vector: activity young (s,ia,is,ih), mid (s,i..), old (s,i..)
objective = function(t, sird, c){
  exp(-(rho+delta))*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[2])
    + sird[4]*u_hosp(c[3]) + sird[5]*u_sus(c[4]) + sird[6]*u_sus(c[4]) + sird[7]*u_sus(c[5])
    + sird[8]*u_hosp(c[6]) + sird[9]*u_sus(c[7]) + sird[10]*u_sus(c[7]) + sird[11]*u_sus(c[8])
    + sird[12]*u_hosp(c[9]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 9)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(0.333, 0.001, 0, 0, 0, 0,
                        0.332, 0.001, 0.001, 0, 0, 0,
                        0.331, 0.001, 0.001, 0.001, 0, 0)),horizon)

sird_init_seed_cond = rep(list(c(0.333, 0.001, 0, 0,
                                 0.332, 0.001, 0.001, 0,
                                 0.331, 0.001, 0.001, 0.001)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(rep(0.8, 9)),horizon))

## used for graphing/charting.
## a: list of activity, a[[t]] is a vector of activities (a_s^y, a_ia^y, ...)_t.
## sird: a pre-initialized-length list of initial time 0
## young (s,ia,is,ih,r,d), mid (s,i..,r,d), old (s,i..,r,d)
soc.trajectory.chart = function(a, sird.init){
  
  sird=sird.init
  
  for(j in 1:(length(a)-1)){
    sird[[j+1]] = next_period_sird(sird[[j]],a[[j]])
  }
  
  sird
}

soc.trajectory = function(a, sird.init){
  
  sird=sird.init
  
  for(j in 1:(length(a)-1)){
    sird[[j+1]] = next_period_sird_cond(sird[[j]],a[[j]])
  }
  
  sird
  
}

soc.value = function(a, sird.init){
  ab = chunk(a, 9)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 500),
                     hessian = FALSE,
                     lower = rep(0.2, horizon),
                     upper = rep(1.1, horizon))
  
  return (x)
  
}
