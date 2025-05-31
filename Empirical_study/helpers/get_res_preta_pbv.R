library(mtrank)
library(magrittr)
library(netmeta)
library(PlackettLuce)
library(tidyverse)

source("Main/helpers/alternative_nma.R")

source("Main/helpers/nma.krahn.output.R")

mod_OR <- readRDS("mod_OR.rds")

mod_OR_new <- list()

id1 <- unique(relax.tcc$id)

mat <- list()

tau <- list()

avg_var <- list()

rel_range <- list()

total_sample <- list()

int <- list()

n_int <- list()

for(i in 1:length(mod_OR)){
  
if(i %in% id1){
  
mod_OR_new[[i]] <- mod_OR[[i]]  

tau[[i]] <- mod_OR_new[[i]]$tau

mat[[i]] <- mod_OR_new[[i]]$seTE.random

## keep all elements except diagonal

mat[[i]] <- mat[[i]][row(mat[[i]]) != col(mat[[i]])]

## transform to variances
mat[[i]] <- mat[[i]]^2

## calculate the average variance
avg_var[[i]] <- mean(mat[[i]])

## calculate the relative range
rel_range[[i]] <- (max(mat[[i]])-min(mat[[i]]))/max(mat[[i]])
                                                    
## total sample size over number of interventions

total_sample[[i]] <- sum(mod_OR_new[[i]]$n1+mod_OR_new[[i]]$n2)

int[[i]] <- length(mod_OR_new[[i]]$trts)

n_int[[i]] <- total_sample[[i]]/int[[i]]

}else{

mod_OR_new[[i]] <- NA

mat[[i]] <- NA

avg_var[[i]] <- NA

rel_range[[i]] <- NA

total_sample[[i]] <- NA

int[[i]] <- NA

n_int[[i]] <- NA

tau[[i]] <- NA

}

  }

mod_OR_new1 <- Filter(function(x) !is.null(x) && !(is.atomic(x) && all(is.na(x))),mod_OR_new)

avg_var1 <- Filter(function(x) !is.null(x) && !(is.atomic(x) && all(is.na(x))),avg_var)

rel_range1 <- Filter(function(x) !is.null(x) && !(is.atomic(x) && all(is.na(x))),rel_range)

n_int1 <- Filter(function(x) !is.null(x) && !(is.atomic(x) && all(is.na(x))),n_int)

tau1 <- Filter(function(x) !is.null(x) && !(is.atomic(x) && all(is.na(x))),tau)

avg_var1 <- unlist(avg_var1)

rel_range1 <- unlist(rel_range1)

n_int1 <- unlist(n_int1)

tau1 <- unlist(tau1)

dat_f <- cbind.data.frame(id1,star1,tau1,cor1,avg_var1,rel_range1,n_int1)

names(dat_f) <- c("id","is.star","tau","corr","avg_var","rel_var_range","n_int")


# write.csv(dat_f,"C:/Users/evrenogl/Desktop/nmadb_ranking/Main/Data/NMA_estimates/nrelax_corrs_metrics.csv",row.names = F)


# 




summary(dat_f$corr)

dat_f1 <- dat_f %>% 
  filter(is.star == "Not a star network")

summary(dat_f1$corr)


saveRDS(mod_OR_new1, file = "mod_OR_relax.rds")

ranks <- list()

small_values_preta <- list()

preta <- list()

res_preta <- list()

for(i in 1:length(mod_OR_new1)){
  
ranks[[i]] <- rankogram(mod_OR_new1[[i]])$ranking.matrix.random[,1]  

small_values_preta[[i]] <- ifelse(mod_OR_new1[[i]]$small.values=="undesirable","bad","good")

preta[[i]] <- alternativenma(mod_OR_new1[[i]],small.values = small_values_preta[[i]])  

res_preta[[i]] <- data.frame("trts"=row.names(preta[[i]]$averages),"preta"=preta[[i]]$averages$Pscoreaverage)

res_preta[[i]] <- res_preta[[i]] %>% 
  arrange(trts)
  
}

ranks1 <- list()

for(i in 1:length(ranks)){
  
ranks1[[i]] <- data.frame("trt"=names(ranks[[i]]),
                          "pBV"=unname(ranks[[i]]))  
  

ranks1[[i]] <- ranks1[[i]] %>% 
  arrange(trt)
  
}

ranks1 <- list_rbind(ranks1)

res_preta1 <- list_rbind(res_preta)

res_pbv_preta <- cbind.data.frame(ranks1$pBV,res_preta1$preta)

names(res_pbv_preta) <- c("pBV","preta")

#  check <- read.csv("Main/Data/NMA_estimates/relax_tcc.csv ")
# 
# all.equal(ranks1$trt,res_preta1$trts,check$id)

# write.csv(res_pbv_preta,"C:/Users/evrenogl/Desktop/nmadb_ranking/Main/Data/NMA_estimates/relax_pbv_preta.csv",
#                               row.names = F)

