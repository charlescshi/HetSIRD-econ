library(ggplot2)
library(dplyr)
library(numDeriv)
library(dfoptim)
library(lubridate)
library(foreach)
library(purrr)
library(optimParallel)
library(extrafont)
library(EpiEstim)
library(Cairo)
CairoWin()

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)

## Parameters

hs = 1
hia = 1
his = 1
hih = 0.8
beta = 1/2.68
gamma_a = 1/6
gamma_s = 1/10
gamma_h = 1/23
sigma = 1
rho = 0.05/365
delta = 0.67/365
tau_a_y = 0.79
tau_s_y = 0.206398
tau_h_y = 0.003602
tau_a_m = 0.5
tau_s_m = 0.43545
tau_h_m = 0.06455
tau_a_o = 0.31
tau_s_o = 0.526
tau_h_o = 0.164
pi_y = 1.315*10^(-4)
pi_m = 0.008415
pi_o = 0.0506
pi_a_y = 0
pi_s_y = 4.102*10^(-4)
pi_h_y = 0.013
pi_a_m = 0
pi_s_m = 5.049*10^(-4)
pi_h_m = 0.060
pi_a_o = 0
pi_s_o = 0.0376
pi_h_o = 0.188
v = 31755
psi = 5.078
kappa_a_y = pi_a_y*v
kappa_s_y = pi_s_y*v
kappa_h_y = pi_h_y*v
kappa_a_m = pi_a_m*v
kappa_s_m = pi_s_m*v
kappa_h_m = pi_h_m*v
kappa_a_o = pi_a_o*v
kappa_s_o = pi_s_o*v
kappa_h_o = pi_h_o*v
minimal_enforcement = 0.1554
chi = 127

sigma_yy = 1.2
sigma_ym = 1
sigma_yo = 0.8
sigma_my = sigma_ym
sigma_mm = 1
sigma_mo = 1
sigma_oy = sigma_yo
sigma_om = sigma_mo
sigma_oo = 1.2

sigma_hy = 0.7
sigma_hm = 0.7
sigma_ho = 0.6

## Initial conditions

uspop = 332599000
initr0 = 3.73

y.d0 = 5/uspop
y.r0 = (y.d0/pi_y) - y.d0
y.i0 = 1 - y.r0 - y.d0 - exp(-initr0*(y.r0 + y.d0))
y.ia0 = y.i0*tau_a_y
y.is0 = y.i0*tau_s_y
y.ih0 = y.i0*tau_h_y
y.s0 = 0.3844 - y.ia0 - y.is0 - y.ih0 - y.r0 - y.d0

m.d0 = 131/uspop
m.r0 = (m.d0/pi_m) - m.d0
m.i0 = 1 - m.r0 - m.d0 - exp(-initr0*(m.r0 + m.d0))
m.ia0 = m.i0*tau_a_m
m.is0 = m.i0*tau_s_m
m.ih0 = m.i0*tau_h_m
m.s0 = 0.4524 - m.ia0 - m.is0 - m.ih0 - m.r0 - m.d0

o.d0 = 440/uspop
o.r0 = (o.d0/pi_o) - o.d0
o.i0 = 1 - o.r0 - o.d0 - exp(-initr0*(o.r0 + o.d0))
o.ia0 = o.i0*tau_a_o
o.is0 = o.i0*tau_s_o
o.ih0 = o.i0*tau_h_o
o.s0 = 1 - o.ia0 - o.is0 - o.ih0 - o.r0 - o.d0 - 0.3844 - 0.4524

## Global settings

horizon = 1000
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
  exp(-(rho+delta)*t)*(
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

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
    )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
      sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
    )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[2])
    + sird[4]*u_hosp(c[3]) + sird[5]*u_sus(c[4]) + sird[6]*u_sus(c[4]) + sird[7]*u_sus(c[5])
    + sird[8]*u_hosp(c[6]) + sird[9]*u_sus(c[7]) + sird[10]*u_sus(c[7]) + sird[11]*u_sus(c[8])
    + sird[12]*u_hosp(c[9]))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 9)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.6, 0.1554, 0.1554, 0.65,
                                      0.1554, 0.1554, 0.7, 0.1554, 0.1554)),horizon))

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
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e10),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*9),
                     upper = rep(1, horizon*9))
  
  return (x)
  
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[2]
    vect[4] = pa[3]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[4]
    vect[8] = pa[4]
    vect[9] = pa[5]
    vect[10] = pa[6]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[7]
    vect[14] = pa[7]
    vect[15] = pa[8]
    vect[16] = pa[9]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

##
## Section
## 3.2.3
## 

eval1 <- soc.optim(sird_init_seed_cond)
eval1.par <- chunk(eval1$par, 9)
eval1.par <- insert.activity.dummies(eval1.par)
eval1.traj <-soc.trajectory.chart(eval1.par, sird_init_seed)
eval1.par <- as.data.frame(do.call(rbind, eval1.par))
eval1.traj <- as.data.frame(do.call(rbind, eval1.traj))
write.csv(eval1.par, "eval1par.csv")
write.csv(eval1.traj, "eval1traj.csv")

## PLOTTING

data.model1 <- data.frame(time_vec, eval1.par, eval1.traj)
data.model1 <- data.model1 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.model1$repno <- esti.r1(data.model1)
data.model1c <- head(data.model1, -400)

theme_set(theme_bw())

##soc activity
p1 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Asymptomatic Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for symptomatic
       and hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p1, filename='p1.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p2 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p2, filename='p2.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p3 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p3, filename='p3.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p4 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p4, filename='p4.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p5 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p5, filename='p5.png', dpi=300, type='cairo', width=10, height=5, units='in')

##
## Section
## 3.2.1
## 

eval2.par <- a_no_behav_cond
eval2.par <- insert.activity.dummies(eval2.par)
eval2.traj <-soc.trajectory.chart(eval2.par, sird_init_seed)
eval2.par <- as.data.frame(do.call(rbind, eval2.par))
eval2.traj <- as.data.frame(do.call(rbind, eval2.traj))
write.csv(eval2.par, "eval2par.csv")
write.csv(eval2.traj, "eval2traj.csv")

## PLOTTING

data.model2 <- data.frame(time_vec, eval2.par, eval2.traj)
data.model2 <- data.model2 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

data.model2$repno <- esti.r1(data.model2)
data.model2c <- head(data.model2, -800)

##soc activity
p6 <- ggplot(data.model2c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Asymptomatic Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for symptomatic
       and hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p6, filename='p6.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p7 <- ggplot(data.model2c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p7, filename='p7.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p8 <- ggplot(data.model2c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p8, filename='p8.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p9 <- ggplot(data.model2c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p9, filename='p9.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p10 <- ggplot(data.model2c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p10, filename='p10.png', dpi=300, type='cairo', width=10, height=5, units='in')

##
## Section
## 3.2.4
## 

next_period_sird_cond = function(sird, c){
  next_sird = vector(mode="numeric", length=12)
  temp_young = sigma*(contact(c[1], c[1])*sird[2]
                      + contact(c[1], c[1])*sird[3]
                      + contact(c[1], c[2])*sird[4]
                      + contact(c[1], c[3])*sird[6]
                      + contact(c[1], c[3])*sird[7]
                      + contact(c[1], c[4])*sird[8]
                      + contact(c[1], c[5])*sird[10]
                      + contact(c[1], c[5])*sird[11]
                      + contact(c[1], c[6])*sird[12])
  temp_mid = sigma*(contact(c[3], c[1])*sird[2]
                    + contact(c[3], c[1])*sird[3]
                    + contact(c[3], c[2])*sird[4]
                    + contact(c[3], c[3])*sird[6]
                    + contact(c[3], c[3])*sird[7]
                    + contact(c[3], c[4])*sird[8]
                    + contact(c[3], c[5])*sird[10]
                    + contact(c[3], c[5])*sird[11]
                    + contact(c[3], c[6])*sird[12])
  temp_old = sigma*(contact(c[5], c[1])*sird[2]
                    + contact(c[5], c[1])*sird[3]
                    + contact(c[5], c[2])*sird[4]
                    + contact(c[5], c[3])*sird[6]
                    + contact(c[5], c[3])*sird[7]
                    + contact(c[5], c[4])*sird[8]
                    + contact(c[5], c[5])*sird[10]
                    + contact(c[5], c[5])*sird[11]
                    + contact(c[5], c[6])*sird[12])
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

objective = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[2]) + sird[5]*u_sus(c[3]) + sird[6]*u_sus(c[3]) + sird[7]*u_sus(c[3])
    + sird[8]*u_hosp(c[4]) + sird[9]*u_sus(c[5]) + sird[10]*u_sus(c[5]) + sird[11]*u_sus(c[5])
    + sird[12]*u_hosp(c[6]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o)
    )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
      sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
    )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[2]) + sird[5]*u_sus(c[3]) + sird[6]*u_sus(c[3]) + sird[7]*u_sus(c[3])
    + sird[8]*u_hosp(c[4]) + sird[9]*u_sus(c[5]) + sird[10]*u_sus(c[5]) + sird[11]*u_sus(c[5])
    + sird[12]*u_hosp(c[6]))
}

a_init_guess = flatten_dbl(rep(list(
  c(0.7, minimal_enforcement, 0.7, minimal_enforcement, 0.7, minimal_enforcement)),horizon))

soc.value = function(a, sird.init){
  ab = chunk(a, 6)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[1]
    vect[4] = pa[2]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[3]
    vect[8] = pa[3]
    vect[9] = pa[3]
    vect[10] = pa[4]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[5]
    vect[14] = pa[5]
    vect[15] = pa[5]
    vect[16] = pa[6]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*6),
                     upper = rep(1, horizon*6))
  
  return (x)
  
}

stopCluster(cl)

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)
clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

eval3 <- soc.optim(sird_init_seed_cond)
eval3.par <- chunk(eval3$par, 6)
eval3.par <- insert.activity.dummies(eval3.par)
eval3.traj <-soc.trajectory.chart(eval3.par, sird_init_seed)
eval3.par <- as.data.frame(do.call(rbind, eval3.par))
eval3.traj <- as.data.frame(do.call(rbind, eval3.traj))
write.csv(eval3.par, "eval3par.csv")
write.csv(eval3.traj, "eval3traj.csv")

## PLOTTING

data.model3 <- data.frame(time_vec, eval3.par, eval3.traj)
data.model3 <- data.model3 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

data.model3$repno <- esti.r1(data.model3)
data.model3c <- head(data.model3, -400)

##soc activity
p11 <- ggplot(data.model3c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Non-hospitalized Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for
       hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))

ggsave(p11, filename='p11.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p12 <- ggplot(data.model3c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')

ggsave(p12, filename='p12.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p13 <- ggplot(data.model3c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p13, filename='p13.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p14 <- ggplot(data.model3c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p14, filename='p14.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p15 <- ggplot(data.model3c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))

ggsave(p15, filename='p15.png', dpi=300, type='cairo', width=10, height=5, units='in')

##
## Section
## 3.2.5
## 

next_period_sird_cond = function(sird, c){
  next_sird = vector(mode="numeric", length=12)
  temp_young = sigma*(contact(c[1], c[1])*sird[2]
                      + contact(c[1], c[1])*sird[3]
                      + contact(c[1], c[1])*sird[4]
                      + contact(c[1], c[2])*sird[6]
                      + contact(c[1], c[2])*sird[7]
                      + contact(c[1], c[2])*sird[8]
                      + contact(c[1], c[3])*sird[10]
                      + contact(c[1], c[3])*sird[11]
                      + contact(c[1], c[3])*sird[12])
  temp_mid = sigma*(contact(c[2], c[1])*sird[2]
                    + contact(c[2], c[1])*sird[3]
                    + contact(c[2], c[1])*sird[4]
                    + contact(c[2], c[2])*sird[6]
                    + contact(c[2], c[2])*sird[7]
                    + contact(c[2], c[2])*sird[8]
                    + contact(c[2], c[3])*sird[10]
                    + contact(c[2], c[3])*sird[11]
                    + contact(c[2], c[3])*sird[12])
  temp_old = sigma*(contact(c[3], c[1])*sird[2]
                    + contact(c[3], c[1])*sird[3]
                    + contact(c[3], c[1])*sird[4]
                    + contact(c[3], c[2])*sird[6]
                    + contact(c[3], c[2])*sird[7]
                    + contact(c[3], c[2])*sird[8]
                    + contact(c[3], c[3])*sird[10]
                    + contact(c[3], c[3])*sird[11]
                    + contact(c[3], c[3])*sird[12])
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

objective = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[1]) + sird[5]*u_sus(c[2]) + sird[6]*u_sus(c[2]) + sird[7]*u_sus(c[2])
    + sird[8]*u_hosp(c[2]) + sird[9]*u_sus(c[3]) + sird[10]*u_sus(c[3]) + sird[11]*u_sus(c[3])
    + sird[12]*u_hosp(c[3]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[1]) + sird[5]*u_sus(c[2]) + sird[6]*u_sus(c[2]) + sird[7]*u_sus(c[2])
    + sird[8]*u_hosp(c[2]) + sird[9]*u_sus(c[3]) + sird[10]*u_sus(c[3]) + sird[11]*u_sus(c[3])
    + sird[12]*u_hosp(c[3])
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
    )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
      sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
    )
  )
}

a_init_guess = flatten_dbl(rep(list(c(0.5, 0.5, 0.5)),horizon))

soc.value = function(a, sird.init){
  ab = chunk(a, 3)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[1]
    vect[4] = pa[1]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[2]
    vect[8] = pa[2]
    vect[9] = pa[2]
    vect[10] = pa[2]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[3]
    vect[14] = pa[3]
    vect[15] = pa[3]
    vect[16] = pa[3]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e8),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*6),
                     upper = rep(1, horizon*6))
  
  return (x)
  
}

stopCluster(cl)

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)
clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

eval4 <- soc.optim(sird_init_seed_cond)
eval4.par <- chunk(eval4$par, 3)
eval4.par <- insert.activity.dummies(eval4.par)
eval4.traj <-soc.trajectory.chart(eval4.par, sird_init_seed)
eval4.par <- as.data.frame(do.call(rbind, eval4.par))
eval4.traj <- as.data.frame(do.call(rbind, eval4.traj))
write.csv(eval4.par, "eval4par.csv")
write.csv(eval4.traj, "eval4traj.csv")

## PLOTTING

data.model4 <- data.frame(time_vec, eval4.par, eval4.traj)
data.model4 <- data.model4 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

data.model4$repno <- esti.r1(data.model4)
data.model4c <- head(data.model4, -400)

##soc activity
p16 <- ggplot(data.model4c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Infected Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))

ggsave(p16, filename='p16.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p17 <- ggplot(data.model4c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')

ggsave(p17, filename='p17.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p18 <- ggplot(data.model4c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p18, filename='p18.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p19 <- ggplot(data.model4c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p19, filename='p19.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p20 <- ggplot(data.model4c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))

ggsave(p20, filename='p20.png', dpi=300, type='cairo', width=10, height=5, units='in')




## -------------------------------------------------------
##
## Section 3.5.2: chi = 1/123
##
## -------------------------------------------------------

## Auxillary function

chi = 5/123

## Disease dynamics

## sird_vector: young (s,ia,is,ih,r,d), mid (s,i..,r,d), old (s,i..,r,d)
## contacts_vector: activity young (s,ia,is,ih,r,d), mid (s,i..,r,d), old (s,i..,r,d)
## this version should be used to chart out the trajectory.
next_period_sird = function(sird, c){
  next_sird = vector(mode="numeric", length=18)
  temp_young = sigma*((1-c[19])*contact(c[1], c[1])*(sird[2] + sird[3])
                      + c[19]*contact(c[1], c[2])*sird[2]
                      + c[19]*contact(c[1], c[3])*sird[3]
                      + contact(c[1], c[4])*sird[4]
                      + (1-c[20])*contact(c[1], c[7])*(sird[8] + sird[9])
                      + c[20]*contact(c[1], c[8])*sird[8]
                      + c[20]*contact(c[1], c[9])*sird[9]
                      + contact(c[1], c[10])*sird[10]
                      + (1-c[21])*contact(c[1], c[13])*(sird[14] + sird[15])
                      + c[21]*contact(c[1], c[14])*sird[14]
                      + c[21]*contact(c[1], c[15])*sird[15]
                      + contact(c[1], c[16])*sird[16])
  temp_mid = sigma*((1-c[19])*contact(c[7], c[1])*(sird[2] + sird[3])
                    + c[19]*contact(c[7], c[2])*sird[2]
                    + c[19]*contact(c[7], c[3])*sird[3]
                    + contact(c[7], c[4])*sird[4]
                    + (1-c[20])*contact(c[7], c[7])*(sird[8] + sird[9])
                    + c[20]*contact(c[7], c[8])*sird[8]
                    + c[20]*contact(c[7], c[9])*sird[9]
                    + contact(c[7], c[10])*sird[10]
                    + (1-c[21])*contact(c[7], c[13])*(sird[14] + sird[15])
                    + c[21]*contact(c[7], c[14])*sird[14]
                    + c[21]*contact(c[7], c[15])*sird[15]
                    + contact(c[7], c[16])*sird[16])
  temp_old = sigma*((1-c[19])*contact(c[13], c[1])*(sird[2] + sird[3])
                    + c[19]*contact(c[13], c[2])*sird[2]
                    + c[19]*contact(c[13], c[3])*sird[3]
                    + contact(c[13], c[4])*sird[4]
                    + (1-c[20])*contact(c[13], c[7])*(sird[8] + sird[9])
                    + c[20]*contact(c[13], c[8])*sird[8]
                    + c[20]*contact(c[13], c[9])*sird[9]
                    + contact(c[13], c[10])*sird[10]
                    + (1-c[21])*contact(c[13], c[13])*(sird[14] + sird[15])
                    + c[21]*contact(c[13], c[14])*sird[14]
                    + c[21]*contact(c[13], c[15])*sird[15]
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
  temp_young = sigma*((1-c[10])*contact(c[1], c[1])*(sird[2] + sird[3])
                      + c[10]*contact(c[1], c[2])*(sird[2] + sird[3])
                      + contact(c[1], c[3])*sird[4]
                      + (1-c[11])*contact(c[1], c[4])*(sird[6] + sird[7])
                      + c[11]*contact(c[1], c[5])*(sird[6] + sird[7])
                      + contact(c[1], c[6])*sird[8]
                      + (1-c[12])*contact(c[1], c[7])*(sird[10] + sird[11])
                      + c[12]*contact(c[1], c[8])*(sird[10] + sird[11])
                      + contact(c[1], c[9])*sird[12])
  temp_mid = sigma*((1-c[10])*contact(c[4], c[1])*(sird[2] + sird[3])
                    + c[10]*contact(c[4], c[2])*(sird[2] + sird[3])
                    + contact(c[4], c[3])*sird[4]
                    + (1-c[11])*contact(c[4], c[4])*(sird[6] + sird[7])
                    + c[11]*contact(c[4], c[5])*(sird[6] + sird[7])
                    + contact(c[4], c[6])*sird[8]
                    + (1-c[12])*contact(c[4], c[7])*(sird[10] + sird[11])
                    + c[12]*contact(c[4], c[8])*(sird[10] + sird[11])
                    + contact(c[4], c[9])*sird[12])
  temp_old = sigma*((1-c[10])*contact(c[7], c[1])*(sird[2] + sird[3])
                    + c[10]*contact(c[7], c[2])*(sird[2] + sird[3])
                    + contact(c[7], c[3])*sird[4]
                    + (1-c[11])*contact(c[7], c[4])*(sird[6] + sird[7])
                    + c[11]*contact(c[7], c[5])*(sird[6] + sird[7])
                    + contact(c[7], c[6])*sird[8]
                    + (1-c[12])*contact(c[7], c[7])*(sird[10] + sird[11])
                    + c[12]*contact(c[7], c[8])*(sird[10] + sird[11])
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
  
  exp(-(rho+delta)*t)*(
    (sird[1] + (1-c[10])*(sird[2]+sird[3]))*u_sus(c[1])
    + c[10]*(sird[2]+sird[3])*u_sus(c[2])
    + sird[4]*u_hosp(c[3])
    + (sird[5] + (1-c[11])*(sird[6]+sird[7]))*u_sus(c[4])
    + c[11]*(sird[6]+sird[7])*u_sus(c[5])
    + sird[8]*u_hosp(c[6])
    + (sird[9] + (1-c[12])*(sird[10]+sird[11]))*u_sus(c[7])
    + c[12]*(sird[10]+sird[11])*u_sus(c[8])
    + sird[12]*u_hosp(c[9])
    - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    ) - chi*(c[10]*(sird[2]+sird[3]+sird[1])
             +c[11]*(sird[5]+sird[6]+sird[7])
             +c[12]*(sird[9]+sird[10]+sird[11]))
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  
  exp(-(rho+delta)*t)*((sird[1] + (1-c[10])*(sird[2]+sird[3]))*u_sus(c[1])
  + c[10]*(sird[2]+sird[3])*u_sus(c[2])
  + sird[4]*u_hosp(c[3])
  + (sird[5] + (1-c[11])*(sird[6]+sird[7]))*u_sus(c[4])
  + c[11]*(sird[6]+sird[7])*u_sus(c[5])
  + sird[8]*u_hosp(c[6])
  + (sird[9] + (1-c[12])*(sird[10]+sird[11]))*u_sus(c[7])
  + c[12]*(sird[10]+sird[11])*u_sus(c[8])
  + sird[12]*u_hosp(c[9]))
}

objective_testing_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(-chi*(c[10]*(sird[2]+sird[3]+sird[1])
                             +c[11]*(sird[5]+sird[6]+sird[7])
                             +c[12]*(sird[9]+sird[10]+sird[11])))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(c(rep(1, 18), 0, 0, 0)),horizon)

a_no_behav_cond = rep(list(c(rep(1, 9), 0, 0, 0)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)


reduction = 5

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.6, 0.1554, 0.1554, 0.65,
                                      0.1554, 0.1554, 0.7, 0.1554, 0.1554, 0.3, 0.4, 0.2)),horizon))

a_init_guess_reduc = flatten_dbl(rep(list(c(1, 0.1554, 0.1554, 1,
                                      0.1554, 0.1554, 1, 0.1554, 0.1554, 0.3, 1, 1)),horizon/reduction))

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


soc.trajectory.reduc = function(a, sird.init){
  sird=sird.init
  
  for(j in 1:(length(a)-1)){
    sird[[5*j-3]] = next_period_sird_cond(sird[[5*j-4]],a[[j]])
    sird[[5*j-2]] = next_period_sird_cond(sird[[5*j-3]],a[[j]])
    sird[[5*j-1]] = next_period_sird_cond(sird[[5*j-2]],a[[j]])
    sird[[5*j]] = next_period_sird_cond(sird[[5*j-1]],a[[j]])
    sird[[5*j+1]] = next_period_sird_cond(sird[[5*j]],a[[j]])
  }
  
  sird
}

soc.value = function(a, sird.init){
  ab = chunk(a, 12)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

soc.value.reduc = function(a, sird.init){
  ab = chunk(a, 12)
  penalty = 1
  for(j in 1:(length(ab)-1)){
    if(j>1){
      penalty = penalty + (((ab[[j]][10] - ab[[j-1]][10])^2
                           + (ab[[j]][11] - ab[[j-1]][11])^2
                           + (ab[[j]][12] - ab[[j-1]][12])^2))/100
    }
  }
  s = soc.trajectory.reduc(ab, sird.init)
  ab = rep(ab, each=reduction)
  ret = Reduce("+",pmap(list(time_vec, s, ab), objective))
  ret/penalty
}





soc.optim.reduc <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess_reduc,
                     fn = soc.value.reduc,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(c(rep(minimal_enforcement, 9), 0, 0, 0), horizon/reduction),
                     upper = rep(c(rep(1, 9), 1, 1, 1), horizon/reduction))
  
  return (x)
  
}

insert.activity.dummies <- function(pars, sird){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    srd = sird[[j]]
    vect = vector(mode="numeric", length=21)
    vect[1] = pa[1]
    vect[2] = pa[2]
    vect[3] = pa[2]
    vect[4] = pa[3]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[4]
    vect[8] = pa[5]
    vect[9] = pa[5]
    vect[10] = pa[6]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[7]
    vect[14] = pa[8]
    vect[15] = pa[8]
    vect[16] = pa[9]
    vect[17] = 1
    vect[18] = 1
    vect[19] = pa[10]
    vect[20] = pa[11]
    vect[21] = pa[12]
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

stopCluster(cl)

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)
clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

init5 <- soc.optim.reduc(sird_init_seed_cond)
initguess5 <- flatten_dbl(rep(chunk(init5$par, 12), each=reduction))

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = initguess5,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(c(rep(minimal_enforcement, 9), 0, 0, 0), horizon),
                     upper = rep(c(rep(1, 9), 1, 1, 1), horizon))
  
  return (x)
  
}

stopCluster(cl)

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)
clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

eval5 <- soc.optim(sird_init_seed_cond)
eval5.traj <- soc.trajectory(chunk(eval5$par, 12), sird_init_seed_cond)
eval5.par <- chunk(eval5$par, 12)
eval5.par <- insert.activity.dummies(eval5.par, eval5.traj)
eval5.traj <-soc.trajectory.chart(eval5.par, sird_init_seed)
eval5.par <- as.data.frame(do.call(rbind, eval5.par))
eval5.traj <- as.data.frame(do.call(rbind, eval5.traj))
write.csv(eval5.par, "eval5par.csv")
write.csv(eval5.traj, "eval5traj.csv")

## PLOTTING

data.model5 <- data.frame(time_vec, eval5.par, eval5.traj)
data.model5 <- data.model5 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1) %>%
  mutate(ty = V19) %>%
  mutate(tm = V20) %>%
  mutate(to = V21)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.model5$repno <- esti.r1(data.model5)
data.model5c <- head(data.model5, -400)

theme_set(theme_bw())

##soc activity
p21 <- ggplot(data.model5c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1)+
  labs(title = "Social Activity of Susceptible and Identified Individuals",
       caption = "Note: The social planner chooses the minimum social activity level for identified individuals.",
       y = "Social activity",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p21, filename='p21.png', dpi=300, type='cairo', width=10, height=5, units='in')

p21a <- ggplot(data.model5c, aes(x=time_vec)) +
  geom_line(aes(y=V19, col = "Young"), size = 1) +
  geom_line(aes(y=V20, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V21, col = "Old"), size = 1)+
  labs(title = "Test Rate for individuals",
       y = "Test Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p21a, filename='p21a.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p22 <- ggplot(data.model5c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p22, filename='p22.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p23 <- ggplot(data.model5c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p23, filename='p23.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p24 <- ggplot(data.model5c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p24, filename='p24.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p25 <- ggplot(data.model5c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p25, filename='p25.png', dpi=300, type='cairo', width=10, height=5, units='in')

## -------------------------------------------------------
##
## Section 3.6: homophily
##
## -------------------------------------------------------

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
  temp_young = sigma_yy*(contact(c[1], c[2])*sird[2]
                      + contact(c[1], c[3])*sird[3]
                      + contact(c[1], c[4])*sird[4])
                      + sigma_ym*(contact(c[1], c[8])*sird[8]
                      + contact(c[1], c[9])*sird[9]
                      + contact(c[1], c[10])*sird[10])
                      + sigma_yo*(contact(c[1], c[14])*sird[14]
                      + contact(c[1], c[15])*sird[15]
                      + contact(c[1], c[16])*sird[16])
  temp_mid = sigma_my*(contact(c[7], c[2])*sird[2]
                    + contact(c[7], c[3])*sird[3]
                    + contact(c[7], c[4])*sird[4])
                    + sigma_mm*(contact(c[7], c[8])*sird[8]
                    + contact(c[7], c[9])*sird[9]
                    + contact(c[7], c[10])*sird[10])
                    + sigma_mo*(contact(c[7], c[14])*sird[14]
                    + contact(c[7], c[15])*sird[15]
                    + contact(c[7], c[16])*sird[16])
  temp_old = sigma_oy*(contact(c[13], c[2])*sird[2]
                    + contact(c[13], c[3])*sird[3]
                    + contact(c[13], c[4])*sird[4])
                    + sigma_om*(contact(c[13], c[8])*sird[8]
                    + contact(c[13], c[9])*sird[9]
                    + contact(c[13], c[10])*sird[10])
                    + sigma_oo*(contact(c[13], c[14])*sird[14]
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
  temp_young = sigma_yy*(contact(c[1], c[1])*sird[2]
                      + contact(c[1], c[2])*sird[3]
                      + contact(c[1], c[3])*sird[4])
                      + sigma_ym*(contact(c[1], c[4])*sird[6]
                      + contact(c[1], c[5])*sird[7]
                      + contact(c[1], c[6])*sird[8])
                      + sigma_yo*(contact(c[1], c[7])*sird[10]
                      + contact(c[1], c[8])*sird[11]
                      + contact(c[1], c[9])*sird[12])
  temp_mid = sigma_my*(contact(c[4], c[1])*sird[2]
                    + contact(c[4], c[2])*sird[3]
                    + contact(c[4], c[3])*sird[4])
                    + sigma_mm*(contact(c[4], c[4])*sird[6]
                    + contact(c[4], c[5])*sird[7]
                    + contact(c[4], c[6])*sird[8])
                    + sigma_mo*(contact(c[4], c[7])*sird[10]
                    + contact(c[4], c[8])*sird[11]
                    + contact(c[4], c[9])*sird[12])
  temp_old = sigma_oy*(contact(c[7], c[1])*sird[2]
                    + contact(c[7], c[2])*sird[3]
                    + contact(c[7], c[3])*sird[4])
                    + sigma_om*(contact(c[7], c[4])*sird[6]
                    + contact(c[7], c[5])*sird[7]
                    + contact(c[7], c[6])*sird[8])
                    + sigma_oo*(contact(c[7], c[7])*sird[10]
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
  exp(-(rho+delta)*t)*(
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

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[2])
    + sird[4]*u_hosp(c[3]) + sird[5]*u_sus(c[4]) + sird[6]*u_sus(c[4]) + sird[7]*u_sus(c[5])
    + sird[8]*u_hosp(c[6]) + sird[9]*u_sus(c[7]) + sird[10]*u_sus(c[7]) + sird[11]*u_sus(c[8])
    + sird[12]*u_hosp(c[9]))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 9)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.6, 0.1554, 0.1554, 0.65,
                                      0.1554, 0.1554, 0.7, 0.1554, 0.1554)),horizon))

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
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e10),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*9),
                     upper = rep(1, horizon*9))
  
  return (x)
  
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[2]
    vect[4] = pa[3]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[4]
    vect[8] = pa[4]
    vect[9] = pa[5]
    vect[10] = pa[6]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[7]
    vect[14] = pa[7]
    vect[15] = pa[8]
    vect[16] = pa[9]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

##
## Section
## 3.6
## 

eval7 <- soc.optim(sird_init_seed_cond)
eval7.par <- chunk(eval7$par, 9)
eval7.par <- insert.activity.dummies(eval7.par)
eval7.traj <-soc.trajectory.chart(eval7.par, sird_init_seed)
eval7.par <- as.data.frame(do.call(rbind, eval7.par))
eval7.traj <- as.data.frame(do.call(rbind, eval7.traj))
write.csv(eval7.par, "eval7par.csv")
write.csv(eval7.traj, "eval7traj.csv")

## PLOTTING

data.model7 <- data.frame(time_vec, eval7.par, eval7.traj)
data.model7 <- data.model7 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.model7$repno <- esti.r1(data.model7)
data.model7c <- head(data.model7, -400)

theme_set(theme_bw())

##soc activity
p31 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Asymptomatic Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for symptomatic
       and hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p31, filename='p31.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p32 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p32, filename='p32.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p33 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p33, filename='p33.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p34 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p34, filename='p34.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p35 <- ggplot(data.model1c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p35, filename='p35.png', dpi=300, type='cairo', width=10, height=5, units='in')

## homophily + information set 2



next_period_sird_cond = function(sird, c){
  next_sird = vector(mode="numeric", length=12)
  temp_young = sigma_yy*(contact(c[1], c[1])*sird[2]
                      + contact(c[1], c[1])*sird[3]
                      + contact(c[1], c[2])*sird[4])
                      + sigma_ym*(contact(c[1], c[3])*sird[6]
                      + contact(c[1], c[3])*sird[7]
                      + contact(c[1], c[4])*sird[8])
                      + sigma_yo*(contact(c[1], c[5])*sird[10]
                      + contact(c[1], c[5])*sird[11]
                      + contact(c[1], c[6])*sird[12])
  temp_mid = sigma_my*(contact(c[3], c[1])*sird[2]
                    + contact(c[3], c[1])*sird[3]
                    + contact(c[3], c[2])*sird[4])
                    + sigma_mm*(contact(c[3], c[3])*sird[6]
                    + contact(c[3], c[3])*sird[7]
                    + contact(c[3], c[4])*sird[8])
                    + sigma_mo*(contact(c[3], c[5])*sird[10]
                    + contact(c[3], c[5])*sird[11]
                    + contact(c[3], c[6])*sird[12])
  temp_old = sigma_oy*(contact(c[5], c[1])*sird[2]
                    + contact(c[5], c[1])*sird[3]
                    + contact(c[5], c[2])*sird[4])
                    + sigma_om*(contact(c[5], c[3])*sird[6]
                    + contact(c[5], c[3])*sird[7]
                    + contact(c[5], c[4])*sird[8])
                    + sigma_oo*(contact(c[5], c[5])*sird[10]
                    + contact(c[5], c[5])*sird[11]
                    + contact(c[5], c[6])*sird[12])
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

objective = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[2]) + sird[5]*u_sus(c[3]) + sird[6]*u_sus(c[3]) + sird[7]*u_sus(c[3])
    + sird[8]*u_hosp(c[4]) + sird[9]*u_sus(c[5]) + sird[10]*u_sus(c[5]) + sird[11]*u_sus(c[5])
    + sird[12]*u_hosp(c[6]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[2]) + sird[5]*u_sus(c[3]) + sird[6]*u_sus(c[3]) + sird[7]*u_sus(c[3])
    + sird[8]*u_hosp(c[4]) + sird[9]*u_sus(c[5]) + sird[10]*u_sus(c[5]) + sird[11]*u_sus(c[5])
    + sird[12]*u_hosp(c[6]))
}

a_init_guess = flatten_dbl(rep(list(
  c(0.7, minimal_enforcement, 0.7, minimal_enforcement, 0.7, minimal_enforcement)),horizon))

soc.value = function(a, sird.init){
  ab = chunk(a, 6)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[1]
    vect[4] = pa[2]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[3]
    vect[8] = pa[3]
    vect[9] = pa[3]
    vect[10] = pa[4]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[5]
    vect[14] = pa[5]
    vect[15] = pa[5]
    vect[16] = pa[6]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*6),
                     upper = rep(1, horizon*6))
  
  return (x)
  
}

stopCluster(cl)

use.core <- c(F, F, F, F, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T)
affinity.mask <- sum(use.core*2^((1:length(use.core))-1))

cl <- makeCluster(20)
setDefaultCluster(cl=cl)
clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

eval8 <- soc.optim(sird_init_seed_cond)
eval8.par <- chunk(eval8$par, 6)
eval8.par <- insert.activity.dummies(eval8.par)
eval8.traj <-soc.trajectory.chart(eval8.par, sird_init_seed)
eval8.par <- as.data.frame(do.call(rbind, eval8.par))
eval8.traj <- as.data.frame(do.call(rbind, eval8.traj))
write.csv(eval8.par, "eval8par.csv")
write.csv(eval8.traj, "eval8traj.csv")

## PLOTTING

data.model8 <- data.frame(time_vec, eval8.par, eval8.traj)
data.model8 <- data.model8 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

data.model8$repno <- esti.r1(data.model8)
data.model8c <- head(data.model8, -400)

##soc activity
p36 <- ggplot(data.model8c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Non-hospitalized Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for
       hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))

ggsave(p36, filename='p36.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
p37 <- ggplot(data.model8c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')

ggsave(p37, filename='p37.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
p38 <- ggplot(data.model8c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(p38, filename='p38.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
p39 <- ggplot(data.model8c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(p39, filename='p39.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
p40 <- ggplot(data.model8c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))

ggsave(p40, filename='p40.png', dpi=300, type='cairo', width=10, height=5, units='in')

## minimal enforcement
## simply change the minimal enforcement


## -------------------------------------------------------
##
## Laissez-Faire
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
## contacts_vector: activity young (s), mid (s), old (s)
## this version should be used during optimiaztion.
next_period_sird_cond = function(sird, c){
  next_sird = vector(mode="numeric", length=12)
  temp_young = sigma*(contact(c[1], c[1])*sird[2]
                      + contact(c[1], his)*sird[3]
                      + contact(c[1], hih)*sird[4]
                      + contact(c[1], c[2])*sird[6]
                      + contact(c[1], his)*sird[7]
                      + contact(c[1], hih)*sird[8]
                      + contact(c[1], c[3])*sird[10]
                      + contact(c[1], his)*sird[11]
                      + contact(c[1], hih)*sird[12])
  temp_mid = sigma*(contact(c[2], c[1])*sird[2]
                    + contact(c[2], his)*sird[3]
                    + contact(c[2], hih)*sird[4]
                    + contact(c[2], c[2])*sird[6]
                    + contact(c[2], his)*sird[7]
                    + contact(c[2], hih)*sird[8]
                    + contact(c[2], c[3])*sird[10]
                    + contact(c[2], his)*sird[11]
                    + contact(c[2], hih)*sird[12])
  temp_old = sigma*(contact(c[3], c[1])*sird[2]
                    + contact(c[3], his)*sird[3]
                    + contact(c[3], hih)*sird[4]
                    + contact(c[3], c[2])*sird[6]
                    + contact(c[3], his)*sird[7]
                    + contact(c[3], hih)*sird[8]
                    + contact(c[3], c[3])*sird[10]
                    + contact(c[3], his)*sird[11]
                    + contact(c[3], hih)*sird[12])
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
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[5]*u_sus(c[2]) + sird[6]*u_sus(c[2]) + sird[9]*u_sus(c[3])
    + sird[10]*u_sus(c[3]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective.y = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) - gamma_a*(
      sird[2]*(kappa_a_y)
    ) - gamma_s*(
      sird[3]*(kappa_s_y) 
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi)
    )
  )
}

objective.m = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[5]*u_sus(c[2]) + sird[6]*u_sus(c[2]) - gamma_a*(
    sird[6]*kappa_a_m 
    ) - gamma_s*(
    sird[7]*kappa_s_m
    ) - gamma_h*(
    sird[8]*(kappa_h_m+psi)
    )
  )
}

objective.o = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[9]*u_sus(c[3])
    + sird[10]*u_sus(c[3]) - gamma_a*(sird[10]*kappa_a_o
    ) - gamma_s*(sird[11]*kappa_s_o
    ) - gamma_h*(sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[5]*u_sus(c[2]) + sird[6]*u_sus(c[2]) + sird[9]*u_sus(c[3])
    + sird[10]*u_sus(c[3]))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 3)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.9, 0.5, minimal_enforcement)),horizon))

## current youth guess:

curr_youth_guess = rep(1, horizon)

## current mid guess:

curr_mid_guess = rep(0.6, horizon)

## current old guess:

curr_old_guess = rep(minimal_enforcement, horizon)

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
  ab = chunk(a, 3)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

soc.y.value = function(a, sird.init){
  temp = c(rbind(a, curr_mid_guess, curr_old_guess))
  ab = chunk(temp, 3)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective.y))
}

soc.m.value = function(a, sird.init){
  temp = c(rbind(curr_youth_guess, a, curr_old_guess))
  ab = chunk(temp, 3)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective.m))
}

soc.o.value = function(a, sird.init){
  temp = c(rbind(curr_youth_guess, curr_mid_guess, a))
  ab = chunk(temp, 3)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective.o))
}

soc.optim.y <- function(sirdinit){
  
  x <- optimParallel(par = curr_youth_guess,
                     fn = soc.y.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e10),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon),
                     upper = rep(1, horizon))
  
  return (x)
  
}

soc.optim.m <- function(sirdinit){
  
  x <- optimParallel(par = curr_mid_guess,
                     fn = soc.m.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e10),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon),
                     upper = rep(1, horizon))
  
  return (x)
  
}

soc.optim.o <- function(sirdinit){
  
  x <- optimParallel(par = curr_old_guess,
                     fn = soc.o.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e10),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon),
                     upper = rep(1, horizon))
  
  return (x)
  
}

clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e10),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*3),
                     upper = rep(1, horizon*3))
  
  return (x)
  
}

test1 <- soc.optim.y(sird_init_seed_cond)

curr_youth_guess <- test1$par

clusterExport(cl, varlist=ls(globalenv()))

test2 <- soc.optim.m(sird_init_seed_cond)

curr_mid_guess <- test2$par

clusterExport(cl, varlist=ls(globalenv()))

test3 <- soc.optim.o(sird_init_seed_cond)

curr_old_guess <- test3$par

for(j in 1:500){
  clusterExport(cl, varlist=ls(globalenv()))
  t1 <- soc.optim.y(sird_init_seed_cond)
  curr_youth_guess <- t1$par
  clusterExport(cl, varlist=ls(globalenv()))
  t2 <- soc.optim.m(sird_init_seed_cond)
  curr_mid_guess <- t2$par
  clusterExport(cl, varlist=ls(globalenv()))
  t3 <- soc.optim.o(sird_init_seed_cond)
  curr_old_guess <- t3$par
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = his
    vect[4] = hih
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[2]
    vect[8] = pa[2]
    vect[9] = his
    vect[10] = hih
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[3]
    vect[14] = pa[3]
    vect[15] = his
    vect[16] = hih
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

##
## Section
## 3.2.3
## 

evalfaire <- soc.optim(sird_init_seed_cond)
evalfaire.par <- chunk(evalfaire$par, 3)
evalfaire.par <- insert.activity.dummies(evalfaire.par)
evalfaire.traj <-soc.trajectory.chart(evalfaire.par, sird_init_seed)
evalfaire.par <- as.data.frame(do.call(rbind, evalfaire.par))
evalfaire.traj <- as.data.frame(do.call(rbind, evalfaire.traj))
write.csv(evalfaire.par, "evalfairepar.csv")
write.csv(evalfaire.traj, "evalfairetraj.csv")

## PLOTTING

data.modelfaire <- data.frame(time_vec, evalfaire.par, evalfaire.traj)
data.modelfaire <- data.modelfaire %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.modelfaire$repno <- esti.r1(data.modelfaire)
data.modelfairec <- head(data.modelfaire, -400)

theme_set(theme_bw())

##soc activity
faire1 <- ggplot(data.modelfairec, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible Individuals",
       caption = "Note: Infected individuals have a social activity level of 1.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(faire1, filename='faire1.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
faire2 <- ggplot(data.modelfairec, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(faire2, filename='faire2.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
faire3 <- ggplot(data.modelfairec, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(faire3, filename='faire3.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
faire4 <- ggplot(data.modelfairec, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(faire4, filename='faire4.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
faire5 <- ggplot(data.modelfairec, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(faire5, filename='faire5.png', dpi=300, type='cairo', width=10, height=5, units='in')



## -------------------------------------------------------
##
## Homogenous
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
                      + contact(c[1], c[1])*sird[6]
                      + contact(c[1], c[2])*sird[7]
                      + contact(c[1], c[3])*sird[8]
                      + contact(c[1], c[1])*sird[10]
                      + contact(c[1], c[2])*sird[11]
                      + contact(c[1], c[3])*sird[12])
  temp_mid = sigma*(contact(c[1], c[1])*sird[2]
                    + contact(c[1], c[2])*sird[3]
                    + contact(c[1], c[3])*sird[4]
                    + contact(c[1], c[1])*sird[6]
                    + contact(c[1], c[2])*sird[7]
                    + contact(c[1], c[3])*sird[8]
                    + contact(c[1], c[1])*sird[10]
                    + contact(c[1], c[2])*sird[11]
                    + contact(c[1], c[3])*sird[12])
  temp_old = sigma*(contact(c[1], c[1])*sird[2]
                    + contact(c[1], c[2])*sird[3]
                    + contact(c[1], c[3])*sird[4]
                    + contact(c[1], c[1])*sird[6]
                    + contact(c[1], c[2])*sird[7]
                    + contact(c[1], c[3])*sird[8]
                    + contact(c[1], c[1])*sird[10]
                    + contact(c[1], c[2])*sird[11]
                    + contact(c[1], c[3])*sird[12])
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
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[2])
    + sird[4]*u_hosp(c[3]) + sird[5]*u_sus(c[1]) + sird[6]*u_sus(c[1]) + sird[7]*u_sus(c[2])
    + sird[8]*u_hosp(c[3]) + sird[9]*u_sus(c[1]) + sird[10]*u_sus(c[1]) + sird[11]*u_sus(c[2])
    + sird[12]*u_hosp(c[3]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[2])
    + sird[4]*u_hosp(c[3]) + sird[5]*u_sus(c[1]) + sird[6]*u_sus(c[1]) + sird[7]*u_sus(c[2])
    + sird[8]*u_hosp(c[3]) + sird[9]*u_sus(c[1]) + sird[10]*u_sus(c[1]) + sird[11]*u_sus(c[2])
    + sird[12]*u_hosp(c[3]))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 9)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.5, 0.1554, 0.1554)),horizon))

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
  ab = chunk(a, 3)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*3),
                     upper = rep(1, horizon*3))
  
  return (x)
  
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[2]
    vect[4] = pa[3]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[1]
    vect[8] = pa[1]
    vect[9] = pa[2]
    vect[10] = pa[3]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[1]
    vect[14] = pa[1]
    vect[15] = pa[2]
    vect[16] = pa[3]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

##
## Section
## 3.2.3
## 

hom1 <- soc.optim(sird_init_seed_cond)
hom1.par <- chunk(hom1$par, 3)
hom1.par <- insert.activity.dummies(hom1.par)
hom1.traj <-soc.trajectory.chart(hom1.par, sird_init_seed)
hom1.par <- as.data.frame(do.call(rbind, hom1.par))
hom1.traj <- as.data.frame(do.call(rbind, hom1.traj))
write.csv(hom1.par, "hom1par.csv")
write.csv(hom1.traj, "hom1traj.csv")

## PLOTTING

data.hom1 <- data.frame(time_vec, hom1.par, hom1.traj)
data.hom1 <- data.hom1 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.hom1$repno <- esti.r1(data.hom1)
data.hom1c <- head(data.hom1, -400)

theme_set(theme_bw())

##soc activity
h1 <- ggplot(data.hom1c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Asymptomatic Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for symptomatic
       and hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h1, filename='h1.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
h2 <- ggplot(data.hom1c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(h2, filename='h2.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
h3 <- ggplot(data.hom1c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(h3, filename='h3.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
h4 <- ggplot(data.hom1c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h4, filename='h4.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
h5 <- ggplot(data.hom1c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h5, filename='h5.png', dpi=300, type='cairo', width=10, height=5, units='in')

## -------------------------------------------------------
##
## Homogenous: Information Set 2
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
                      + contact(c[1], c[1])*sird[3]
                      + contact(c[1], c[2])*sird[4]
                      + contact(c[1], c[1])*sird[6]
                      + contact(c[1], c[1])*sird[7]
                      + contact(c[1], c[2])*sird[8]
                      + contact(c[1], c[1])*sird[10]
                      + contact(c[1], c[1])*sird[11]
                      + contact(c[1], c[2])*sird[12])
  temp_mid = sigma*(contact(c[1], c[1])*sird[2]
                    + contact(c[1], c[1])*sird[3]
                    + contact(c[1], c[2])*sird[4]
                    + contact(c[1], c[1])*sird[6]
                    + contact(c[1], c[1])*sird[7]
                    + contact(c[1], c[2])*sird[8]
                    + contact(c[1], c[1])*sird[10]
                    + contact(c[1], c[1])*sird[11]
                    + contact(c[1], c[2])*sird[12])
  temp_old = sigma*(contact(c[1], c[1])*sird[2]
                    + contact(c[1], c[1])*sird[3]
                    + contact(c[1], c[2])*sird[4]
                    + contact(c[1], c[1])*sird[6]
                    + contact(c[1], c[1])*sird[7]
                    + contact(c[1], c[2])*sird[8]
                    + contact(c[1], c[1])*sird[10]
                    + contact(c[1], c[1])*sird[11]
                    + contact(c[1], c[2])*sird[12])
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
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[2]) + sird[5]*u_sus(c[1]) + sird[6]*u_sus(c[1]) + sird[7]*u_sus(c[1])
    + sird[8]*u_hosp(c[2]) + sird[9]*u_sus(c[1]) + sird[10]*u_sus(c[1]) + sird[11]*u_sus(c[1])
    + sird[12]*u_hosp(c[2]) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c[1]) + sird[2]*u_sus(c[1]) + sird[3]*u_sus(c[1])
    + sird[4]*u_hosp(c[2]) + sird[5]*u_sus(c[1]) + sird[6]*u_sus(c[1]) + sird[7]*u_sus(c[1])
    + sird[8]*u_hosp(c[2]) + sird[9]*u_sus(c[1]) + sird[10]*u_sus(c[1]) + sird[11]*u_sus(c[1])
    + sird[12]*u_hosp(c[2]))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 9)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.5, 0.1554)),horizon))

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
  ab = chunk(a, 2)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*2),
                     upper = rep(1, horizon*2))
  
  return (x)
  
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[1]
    vect[4] = pa[2]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[1]
    vect[8] = pa[1]
    vect[9] = pa[1]
    vect[10] = pa[2]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[1]
    vect[14] = pa[1]
    vect[15] = pa[1]
    vect[16] = pa[2]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

##
## Section
## 3.2.3
## 

hom2 <- soc.optim(sird_init_seed_cond)
hom2.par <- chunk(hom2$par, 2)
hom2.par <- insert.activity.dummies(hom2.par)
hom2.traj <-soc.trajectory.chart(hom2.par, sird_init_seed)
hom2.par <- as.data.frame(do.call(rbind, hom2.par))
hom2.traj <- as.data.frame(do.call(rbind, hom2.traj))
write.csv(hom2.par, "hom2par.csv")
write.csv(hom2.traj, "hom2traj.csv")

## PLOTTING

data.hom2 <- data.frame(time_vec, hom2.par, hom2.traj)
data.hom2 <- data.hom2 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.hom2$repno <- esti.r1(data.hom2)
data.hom2c <- head(data.hom2, -400)

theme_set(theme_bw())

##soc activity
h6 <- ggplot(data.hom2c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Asymptomatic Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for symptomatic
       and hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h6, filename='h6.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
h7 <- ggplot(data.hom2c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(h7, filename='h7.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
h8 <- ggplot(data.hom2c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(h8, filename='h8.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
h9 <- ggplot(data.hom2c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h9, filename='h9.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
h10 <- ggplot(data.hom2c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h10, filename='h10.png', dpi=300, type='cairo', width=10, height=5, units='in')

## -------------------------------------------------------
##
## Homogenous: Information Set 3
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
  temp_young = sigma*(contact(c, c)*sird[2]
                      + contact(c, c)*sird[3]
                      + contact(c, c)*sird[4]
                      + contact(c, c)*sird[6]
                      + contact(c, c)*sird[7]
                      + contact(c, c)*sird[8]
                      + contact(c, c)*sird[10]
                      + contact(c, c)*sird[11]
                      + contact(c, c)*sird[12])
  temp_mid = sigma*(contact(c, c)*sird[2]
                    + contact(c, c)*sird[3]
                    + contact(c, c)*sird[4]
                    + contact(c, c)*sird[6]
                    + contact(c, c)*sird[7]
                    + contact(c, c)*sird[8]
                    + contact(c, c)*sird[10]
                    + contact(c, c)*sird[11]
                    + contact(c, c)*sird[12])
  temp_old = sigma*(contact(c, c)*sird[2]
                    + contact(c, c)*sird[3]
                    + contact(c, c)*sird[4]
                    + contact(c, c)*sird[6]
                    + contact(c, c)*sird[7]
                    + contact(c, c)*sird[8]
                    + contact(c, c)*sird[10]
                    + contact(c, c)*sird[11]
                    + contact(c, c)*sird[12])
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
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c) + sird[2]*u_sus(c) + sird[3]*u_sus(c)
    + sird[4]*u_hosp(c) + sird[5]*u_sus(c) + sird[6]*u_sus(c) + sird[7]*u_sus(c)
    + sird[8]*u_hosp(c) + sird[9]*u_sus(c) + sird[10]*u_sus(c) + sird[11]*u_sus(c)
    + sird[12]*u_hosp(c) - gamma_a*(
      sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
    ) - gamma_s*(
      sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
    ) - gamma_h*(
      sird[4]*(kappa_h_y+psi) + sird[8]*(kappa_h_m+psi) + sird[12]*(kappa_h_o+psi)
    )
  )
}

objective_death_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_a*(
    sird[2]*(kappa_a_y) + sird[6]*kappa_a_m + sird[10]*kappa_a_o
  ) - gamma_s*(
    sird[3]*(kappa_s_y) + sird[7]*kappa_s_m + sird[11]*kappa_s_o
  ) - gamma_h*(
    sird[4]*(kappa_h_y) + sird[8]*(kappa_h_m) + sird[12]*(kappa_h_o)
  )
  )
}

objective_hosp_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(- gamma_h*(
    sird[4]*(psi) + sird[8]*(psi) + sird[12]*(psi)
  )
  )
}

objective_lockdown_costs = function(t, sird, c){
  exp(-(rho+delta)*t)*(
    sird[1]*u_sus(c) + sird[2]*u_sus(c) + sird[3]*u_sus(c)
    + sird[4]*u_hosp(c) + sird[5]*u_sus(c) + sird[6]*u_sus(c) + sird[7]*u_sus(c)
    + sird[8]*u_hosp(c) + sird[9]*u_sus(c) + sird[10]*u_sus(c) + sird[11]*u_sus(c)
    + sird[12]*u_hosp(c))
}

## used for testing w/ no behavioral responses:
a_no_behav = rep(list(rep(1, 18)),horizon)

a_no_behav_cond = rep(list(rep(1, 9)),horizon)

## a random initial sird seed
sird_init_seed = rep(list(c(y.s0, y.ia0, y.is0, y.ih0, y.r0, y.d0,
                            m.s0, m.ia0, m.is0, m.ih0, m.r0, m.d0,
                            o.s0, o.ia0, o.is0, o.ih0, o.r0, o.d0)),horizon)

sird_init_seed_cond = rep(list(c(y.s0, y.ia0, y.is0, y.ih0,
                                 m.s0, m.ia0, m.is0, m.ih0,
                                 o.s0, o.ia0, o.is0, o.ih0)),horizon)

time_vec = 1:horizon

## used for initial guess:
a_init_guess = flatten_dbl(rep(list(c(0.5)),horizon))

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
  ab = chunk(a, 1)
  Reduce("+",pmap(list(time_vec, soc.trajectory(ab, sird.init), ab), objective))
}

clusterExport(cl, varlist=ls(globalenv()))
clusterEvalQ(cl, library("purrr"))

soc.optim <- function(sirdinit){
  
  x <- optimParallel(par = a_init_guess,
                     fn = soc.value,
                     sird.init=sirdinit,
                     method = "BFGS",
                     control = list(fnscale = -1, maxit = 5000, trace = 4,
                                    REPORT = 1, factr=1e9),
                     hessian = FALSE,
                     lower = rep(minimal_enforcement, horizon*1),
                     upper = rep(1, horizon*1))
  
  return (x)
  
}

insert.activity.dummies <- function(pars){
  v = list()
  for(j in 1:length(pars)){
    pa = pars[[j]]
    vect = vector(mode="numeric", length=18)
    vect[1] = pa[1]
    vect[2] = pa[1]
    vect[3] = pa[1]
    vect[4] = pa[1]
    vect[5] = 1
    vect[6] = 1
    vect[7] = pa[1]
    vect[8] = pa[1]
    vect[9] = pa[1]
    vect[10] = pa[1]
    vect[11] = 1
    vect[12] = 1
    vect[13] = pa[1]
    vect[14] = pa[1]
    vect[15] = pa[1]
    vect[16] = pa[1]
    vect[17] = 1
    vect[18] = 1
    v = rbind(v, vect)
  }
  v <- split(v, rep(1:ncol(v), each=nrow(v)))
  v <- transpose(v)
  for(j in 1:length(v)){
    v[[j]] = flatten_dbl(v[[j]])
  }
  v
}

##
## Section
## 3.2.3
## 

hom3 <- soc.optim(sird_init_seed_cond)
hom3.par <- chunk(hom3$par, 1)
hom3.par <- insert.activity.dummies(hom3.par)
hom3.traj <-soc.trajectory.chart(hom3.par, sird_init_seed)
hom3.par <- as.data.frame(do.call(rbind, hom3.par))
hom3.traj <- as.data.frame(do.call(rbind, hom3.traj))
write.csv(hom3.par, "hom3par.csv")
write.csv(hom3.traj, "hom3traj.csv")

## PLOTTING

data.hom3 <- data.frame(time_vec, hom3.par, hom3.traj)
data.hom3 <- data.hom3 %>%
  mutate(infy = (V2.1+V3.1+V4.1)/0.3844) %>%
  mutate(infm = (V8.1+V9.1+V10.1)/0.4524) %>%
  mutate(info = (V14.1+V15.1+V16.1)/0.1632) %>%
  mutate(dy = V6.1/0.3844) %>%
  mutate(dm = V12.1/0.4524) %>%
  mutate(do = V18.1/0.1632) %>%
  mutate(s = V1.1+V7.1+V13.1) %>%
  mutate(i = infy*0.3844+infm*0.4524+info*0.1632) %>%
  mutate(r = V5.1+V11.1+V17.1) %>%
  mutate(d = V6.1+V12.1+V18.1) %>%
  mutate(ia = V2.1+V8.1+V14.1) %>%
  mutate(is = V3.1+V9.1+V15.1) %>%
  mutate(ih = V4.1+V10.1+V16.1)

# estimating r0

esti.r1 = function(datamodel){
  r = vector(mode="numeric", length=nrow(datamodel))
  for(j in 1:nrow(datamodel)){
    l = vector(mode="numeric", length=ceiling(2/gamma_h))
    for(s in 1:length(l)){
      if(j-s<=0){
        l[s] = 0
      }
      else{
        l[s] = datamodel$ia[j-s] * (1 - gamma_a)^s +
          datamodel$is[j-s] * (1 - gamma_s)^s +
          datamodel$ih[j-s] * (1 - gamma_h)^s
      }
    }
    sm = sum(l)
    for(s in 1:length(l)){
      l[s] = l[s]/sm
    }
    z = 0
    for(s in 1:length(l)){
      if(j-s>0){
        z = z + datamodel$i[j-s]*l[s]
      }
    }
    r[j] = datamodel$i[j]/sum(z)
  }
  r[1] = NA
  r
}
data.hom3$repno <- esti.r1(data.hom3)
data.hom3c <- head(data.hom3, -400)

theme_set(theme_bw())

##soc activity
h11 <- ggplot(data.hom3c, aes(x=time_vec)) +
  geom_line(aes(y=V1, col = "Young"), size = 1) +
  geom_line(aes(y=V7, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=V13, col = "Old"), size = 1) +
  labs(title = "Social Activity of Susceptible and Asymptomatic Individuals",
       caption = "Note: The social planner enforces a social activity level of 1 on recovered
       individuals, and enforces the minimum social activity level possible for symptomatic
       and hospitalized individuals.",
       y = "Social activity (optimal level of contacts at 1)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h11, filename='h11.png', dpi=300, type='cairo', width=10, height=5, units='in')

## infections
h12 <- ggplot(data.hom3c, aes(x=time_vec)) +
  geom_line(aes(y=infy, col = "Young"), size = 1) +
  geom_line(aes(y=infm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=info, col = "Old"), size = 1) +
  labs(title = "Infection Rates",
       caption = "Note: Infection Rates are conditional on each subset of population.
       The scaling is log.",
       y = "Infection Rate",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(h12, filename='h12.png', dpi=300, type='cairo', width=10, height=5, units='in')

## deaths
h13 <- ggplot(data.hom3c, aes(x=time_vec)) +
  geom_line(aes(y=dy, col = "Young"), size = 1) +
  geom_line(aes(y=dm, col = "Middle-aged"), size = 1) +
  geom_line(aes(y=do, col = "Old"), size = 1) +
  labs(title = "Cumulative deaths",
       caption = "Note: Cumulative deaths are conditional on each subset of population.
       The scaling is log.",
       y = "Cumulative Deaths",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16)) +
  scale_y_continuous(trans='log10')
ggsave(h13, filename='h13.png', dpi=300, type='cairo', width=10, height=5, units='in')

## overall_sird
h14 <- ggplot(data.hom3c, aes(x=time_vec)) +
  geom_line(aes(y=s, col = "Susceptible"), size = 1) +
  geom_line(aes(y=i, col = "Infected"), size = 1) +
  geom_line(aes(y=r, col = "Recovered"), size = 1) +
  geom_line(aes(y=d, col = "Dead"), size = 1) +
  labs(title = "Disease Trajectory",
       caption = "Aggregate S/I/R/D trajectory for all population groups.",
       y = "Share of Population",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h14, filename='h14.png', dpi=300, type='cairo', width=10, height=5, units='in')

## r(t)
h15 <- ggplot(data.hom3c, aes(x=time_vec)) +
  geom_line(aes(y=repno, col = "R(t)"), size = 1) +
  labs(title = "Effective Reproduction Number",
       caption = "Aggregate Reproduction Number over time.",
       y = "R(t)",
       x = "Time (Days)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.caption = element_text(hjust = 0)) +
  theme(text=element_text(family="Segoe UI", size=16))
ggsave(h15, filename='h15.png', dpi=300, type='cairo', width=10, height=5, units='in')