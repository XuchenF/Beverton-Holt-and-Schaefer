params<-c("h","K", "lnq", "P", "Imed", "ct", "x")
hh=0.66
rr=1.18
KK<-3350000
qq<- -18
jags.inits<-function(){
list("h"=hh, "K"=KK, "lnq"=qq)
}

##Schaefer production function jags model
RJBugs<-function(){
K~dunif(1000000, 5580000)
lnq~dunif(-25, -5)
q<-exp(lnq)
r~dlnorm(0.166, 55)
isigma2 ~ dgamma(2,0.01)
sigma2 <- 1/isigma2+pow(0.15,2)
sigma <- pow(sigma2,0.5)
alpha~dlnorm(log(0.7), (log(0.7)-log(0.55))/4)
for(t in 1:N){
ct[t] ~ dlnorm(log(C[t]),pow(0.15,-2))
}
itau2~dgamma(4,0.01)
eps=0.01
penm[1] <- 0 # no penalty for first biomass
Pmean[1] <- log(alpha)
P[1] ~ dlnorm(Pmean[1],itau2)
for (t in 2:N) {
Pmean[t] <- ifelse(P[t-1] > 0.25,
log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct[t-1]/K,eps)),  # Process equation
log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct[t-1]/K,eps)))
P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*K*P[t])-log(q*K*(eps+0.001)),ifelse(P[t]>1,log(q*K*P[t])-log(q*K*(0.99)),0)) # penalty if Pmean is outside viable biomass
}
for(i in 1:N){
pen.bk[i] ~ dnorm(penm[i],10000)
Imed[i]<-log(q*P[i]*K)
I[i]~dlnorm(Imed[i], pow(sigma2, -1))
}
}


##BHDPF jags model
RJBugs<-function(){
    K~dunif(1000000, 5580000)
    lnq~dunif(-25, -5)
    q<-exp(lnq)
    h~dlnorm(-0.416, 40)
    m=0.59
    isigma2 ~ dgamma(2,0.01)
    sigma2 <- 1/isigma2+pow(0.15,2)
    sigma <- pow(sigma2,0.5)
    alpha~dlnorm(log(0.7), (log(0.7)-log(0.55))/4)
    for(t in 1:N){
        ct[t] ~ dlnorm(log(C[t]),pow(0.15,-2))
    }
    itau2~dgamma(4,0.01)
    eps=0.01
    penm[1] <- 0 # no penalty for first biomass
    Pmean[1] <- log(alpha)
    P[1] ~ dlnorm(Pmean[1],itau2)
    for (t in 2:N) {
        Pmean[t] <- log(max(P[t-1]*(1-m) + h*P[t-1]*4*m/(1+P[t-1]*(h-0.2)/(0.2*K-0.2*h*K)*K) - ct[t-1]/K,eps))
        P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
        penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*K*P[t])-log(q*K*(eps+0.001)),ifelse(P[t]>1,log(q*K*P[t])-log(q*K*(0.99)),0)) # penalty if Pmean is outside viable biomass
    }
    for(i in 1:N){
        pen.bk[i] ~ dnorm(penm[i],10000)
        Imed[i]<-log(q*P[i]*K)
        I[i]~dlnorm(Imed[i], pow(sigma2, -1))
    }
}


##BHDPF2 jags model
RJBugs<-function(){
    K~dunif(1000000, 5580000)
    lnq~dunif(-25, -5)
    q<-exp(lnq)
    h~dlnorm(-0.416, 40)
    m=0.59
    isigma2 ~ dgamma(2,0.01)
    sigma2 <- 1/isigma2+pow(0.15,2)
    sigma <- pow(sigma2,0.5)
    alpha~dlnorm(log(0.7), (log(0.7)-log(0.55))/4)
    for(t in 1:N){
        ct[t] ~ dlnorm(log(C[t]),pow(0.15,-2))
    }
    itau2~dgamma(4,0.01)
    eps=0.01
    penm[1] <- 0 # no penalty for first biomass
    Pmean[1] <- log(alpha)
    P[1] ~ dlnorm(Pmean[1],itau2)
    for (t in 2:N) {
        Pmean[t] <- log(max(P[t-1]*(1-m) + h*P[t-1]*4*m/(1+P[t-1]*(h-0.2)/(0.2*K-0.2*h*K)*K) - ct[t-1]/K,eps))
        P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
        penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*K*P[t])-log(q*K*(eps+0.001)),ifelse(P[t]>1,log(q*K*P[t])-log(q*K*(0.99)),0)) # penalty if Pmean is outside viable biomass
    }
    for(i in 1:N){
        pen.bk[i] ~ dnorm(penm[i],10000)
        Imed[i]<-log(q*P[i]*K)
        I[i]~dlnorm(Imed[i], pow(sigma2, -1))
    }
}
