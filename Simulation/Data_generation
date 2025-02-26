library(tidyverse)
library(netmeta)
library(metafor)
library(dplyr)
library(MASS)

'%!in%' <- function(x,y)!('%in%'(x,y))

source("Simulation/helpers_sim/dlong.R")

#### Settings 
set.seed(444)
N.sim=1000

data1=list()
data2=list()
logOR=list()
logOR1=list()
OR=list()
sample <- list()


NT=8 #### number of treatments in the network
NS=4 #### number of studies per comparison

## Define the heterogeneity in the network

tau=0 ####  heterogeneity SD
#tau=0.1 ####  heterogeneity SD
#tau=0.28 ####  heterogeneity SD
#tau=0.5 ####  heterogeneity SD

#tau=0.7 ####  heterogeneity SD
#tau=1 ####  heterogeneity SD

### Define minimum and maximum number of patiets per treatment arm

Npmin=50 #### minimum number of patients per arm
Npmax=100 #### maximum number of patients per arm

#Npmin=30 #### minimum number of patients per arm
#Npmax=300 #### maximum number of patients per arm

### Define how many of the total generated comparison will be retained.

keep <- 1
#keep <- 0.85
#keep <- 0.75
#keep <- 0.65
#keep <- 0.55

### define treatment indices
t1=c()
t2=c()
for (i in 1:(NT-1)){
  for (k in (i+1):NT){
    for(j in 1:NS){
      t1=c(t1,i)
      t2=c(t2,k)      }}}
N.stud=length(t1)

### define patients per treatment arm
for (i in 1:N.sim)
{   
  logOR[[i]]=seq(from =1/(NT-1), to = 1, by = 1/(NT-1))  ### true logOR across studies 
  OR[[i]]=c(1,exp(logOR[[i]]))
  
  
  data1[[i]]=data.frame(t1,t2)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$n1=data1[[i]]$n2=round(runif(N.stud,Npmin,Npmax))}


#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
  data1[[i]]$p.ref=runif(N.stud,0.2,0.4) #### study-specific probability of an event in treatment 1
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
}

#### define probabilities per treatment, per study arm
Sigma=matrix(c(tau^2,tau^2/2,tau^2/2,tau^2),nrow=2)

for (i in 1:N.sim)
{   
  logOR1[[i]]=c(0,logOR[[i]])
  for(j in 1:N.stud){
    
    
    data1[[i]]$truelogOR.t1[j]=logOR1[[i]][data1[[i]]$t1[j]]
    data1[[i]]$truelogOR.t2[j]=logOR1[[i]][data1[[i]]$t2[j]]
    test1=mvrnorm(1,c(data1[[i]]$truelogOR.t1[j],data1[[i]]$truelogOR.t2[j]),Sigma)
    test2=rnorm(1,data1[[i]]$truelogOR.t2[j],tau)
    data1[[i]]$st.sp.logOR.t1[j]=test1[1]*(data1[[i]]$t1[j]!=1)
    data1[[i]]$st.sp.logOR.t2[j]=test1[2]*(data1[[i]]$t1[j]!=1)+test2*(data1[[i]]$t1[j]==1)
    
    data1[[i]]$odds.t1[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t1[j])
    data1[[i]]$odds.t2[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t2[j])
    data1[[i]]$p.t1[j]=data1[[i]]$odds.t1[j]/(1+data1[[i]]$odds.t1[j])
    data1[[i]]$p.t2[j]=data1[[i]]$odds.t2[j]/(1+data1[[i]]$odds.t2[j])
  }}


#### generate the data
for (i in 1:N.sim)
{  
  for(j in 1:N.stud){
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n1[j],data1[[i]]$p.t1[j]) 
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n2[j],data1[[i]]$p.t2[j]) 
    data1[[i]]$comp[j] = paste(data1[[i]]$t1[j],"-",data1[[i]]$t2[j],sep = "")
  }
  sample[[i]] <- sample(unique(data1[[i]]$comp),size = round(keep*length(unique(data1[[i]]$comp)),digits = 0),replace = F)
  
  data1[[i]] <- data1[[i]] %>%
    filter(comp %in% sample[[i]])
  }

for (i in 1:N.sim){ data1[[i]] <- data1[[i]][,-c(6:12)]
data2[[i]] <- dlong(data1[[i]])
row.names(data2[[i]]) <- NULL
}


count <- list()

rm(data1)

for(i in 1:length(data2)){

if(length(unique(data2[[i]]$treat))!=length(1:NT)){
  
  count[[i]] <- 1

  }else{
  count[[i]] <- 0
}
}

count1 <- unlist(count)

E <- which(count1==1)

length(E)

if(length(E)>0){
  
  data2[E[1:length(E)]] <- NULL  
  
}
