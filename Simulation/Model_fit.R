### Load the necessary libraries, all of which are available in CRAN. 
library(tidyverse)
library(netmeta)
library(metafor)
library(dplyr)
library(MASS)
library(parallel)
library(ircor)
library(mtrank)
library(rlist)

## set a working directory
setwd("./Simulation")

## Read files
source("helpers_sim/average_overlap.R")
source("helpers_sim/nma.krahn.output.R")
source("helpers_sim/alternative_nma.R")


# Initiate cluster
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)

N.sim <- length(data2)

X1=list()

for (i in 1:N.sim){ X1[[i]]=list("data"=data2[[i]],"logOR"=logOR[[i]],"tau"=tau[[1]])}
                    
################################################################
###### Ranking models ###########

models=function(X)
{
  tryCatch({
  
  true_rank <- seq(1:NT)
  
  if((length(true_rank)%% 2) == 0) {
    half_true <- true_rank[1:(NT/2)]
  }else{
    half_true <- true_rank[1:((NT-1)/2)]
    }
  
  ##### Ability based metric vs true
  
  kendall_PL_true <- c()
  spearman_PL_true <- c()
  overlap_PL_true <- c()
  
  ## P-scores (random) vs true
  
  kendall_prandom_true <- c()
  spearman_prandom_true <- c()
  overlap_prandom_true <- c()
  
  ## P-scores (common) vs true
  
  kendall_pcommon_true <- c()
  spearman_pcommon_true <- c()
  overlap_pcommon_true <- c()
  
  ## PRETA vs true
  
  kendall_preta_true <- c()
  spearman_preta_true <- c()
  overlap_preta_true <- c()
  
  ## pBV random vs true
  
  kendall_pBV_random_true <- c()
  spearman_pBV_random_true <- c()
  overlap_pBV_random_true <- c()
   
  ######### MODELS
  
  p <- pairwise(data = X$data,studlab = study,treat = treat,n=n,event = events,sm="OR")
  
  ### mtrank
  
  rankings <- tcc(p,
                  mcid = 1.01,
                  small.values ="desirable" )
                          
  
  mod_mtrank <- mtrank(rankings) 
  
  rank_PL <- mod_mtrank$estimates %>% 
    dplyr::select(treatment,log_ability) %>% 
    arrange(desc(log_ability))
  
  list_PL <- as.numeric(rank_PL$treatment)
  
  #### netmeta 
  
  mod_netmeta <- netmeta(p,reference.group = 1,small.values = "desirable")

  ### P-scores from random effects model 
  A <- netrank(mod_netmeta,random = T,common = F)

  p_random <- as.data.frame(A$ranking.random)

  names(p_random) <- "P_score"

  p_random <- p_random %>% 
    arrange(desc(P_score))

  list_p_random <- as.numeric(rownames(p_random))

  ### P-scores from common effect model 
  B <- netrank(mod_netmeta,random = F,fixed = T)
  
  p_common <- as.data.frame(B$ranking.fixed)
  
  names(p_common) <- "P_score"

  p_common <- p_common %>% 
    arrange(desc(P_score))
  
  list_p_common <- as.numeric(rownames(p_common))
  
  #### PRETA
  
  mod_preta <- alternativenma(mod_netmeta,small.values = "good")
  
  res_preta <- mod_preta$averages
  
  res_preta <- res_preta %>% 
    arrange(desc(Pscoreaverage))
  
  list_preta <- as.numeric(row.names(res_preta))
    
  ##### pBV 
  
  rankogram <- rankogram(mod_netmeta)
    
  ### random effects
  pbV_random <- data.frame("pbV_random" = rankogram$ranking.matrix.random[,1])
  
  pbV_random <- pbV_random %>% 
    arrange(desc(pbV_random))
  
  pbV_random$treat <- row.names(pbV_random)
  
  list_pbv_random <- as.numeric(pbV_random$treat)
  
  
  ######## Performance measures
  
  #### similarity Ability based metrics vs true
  
  kendall_PL_true <- c(kendall_PL_true, (cor(list_PL, true_rank, method = "kendall")))
  
  spearman_PL_true <- c(spearman_PL_true, (cor(list_PL, true_rank, method = "spearman")))
  
  overlap_PL_true <- c(overlap_PL_true,(averageoverlap(x=list_PL, y=true_rank,k=length(half_true))))
  
  #### similarity P-scores random vs true 
  
  kendall_prandom_true <- c(kendall_prandom_true, (cor(list_p_random, true_rank, method = "kendall")))
  
  spearman_prandom_true <- c(spearman_prandom_true, (cor(list_p_random, true_rank, method = "spearman")))
  
  overlap_prandom_true <- c(overlap_prandom_true,(averageoverlap(x=list_p_random, y=true_rank,k=length(half_true))))

  #### similarity P-scores common vs true 
  
  kendall_pcommon_true <- c(kendall_pcommon_true, (cor(list_p_common, true_rank, method = "kendall")))
  
  spearman_pcommon_true <- c(spearman_pcommon_true, (cor(list_p_common, true_rank, method = "spearman")))
  
  overlap_pcommon_true <- c(overlap_pcommon_true,(averageoverlap(x=list_p_common, y=true_rank,k=length(half_true))))
  
  #### similarity PRETA vs true
  
  kendall_preta_true <- c(kendall_preta_true, (cor(list_preta, true_rank, method = "kendall")))
  
  spearman_preta_true <- c(spearman_preta_true, (cor(list_preta, true_rank, method = "spearman")))
  
  overlap_preta_true <- c(overlap_preta_true,(averageoverlap(x=list_preta, y=true_rank,k=length(half_true))))
  
  #### similarity pbv random vs true 
  
  kendall_pBV_random_true <- c(kendall_pBV_random_true, (cor(list_pbv_random, true_rank, method = "kendall")))
  
  spearman_pBV_random_true <- c(spearman_pBV_random_true, (cor(list_pbv_random, true_rank, method = "spearman")))
  
  overlap_pBV_random_true <- c(overlap_pBV_random_true,(averageoverlap(x=list_pbv_random, y=true_rank,k=length(half_true))))
  
  
  return(list(
              "kendall_PL_true" = kendall_PL_true,
              "spearman_PL_true" = spearman_PL_true,
              "overlap_PL_true" = overlap_PL_true,
              "kendall_prandom_true" = kendall_prandom_true,
              "spearman_prandom_true" = spearman_prandom_true,
              "overlap_prandom_true" = overlap_prandom_true,
              "kendall_pcommon_true" = kendall_pcommon_true,
              "spearman_pcommon_true" = spearman_pcommon_true,
              "overlap_pcommon_true" = overlap_pcommon_true,
              "kendall_pBV_random_true " = kendall_pBV_random_true, 
              "spearman_pBV_random_true" = spearman_pBV_random_true,
              "overlap_pBV_random_true" = overlap_pBV_random_true, 
              "kendall_preta_true" = kendall_preta_true,
              "spearman_preta_true" = spearman_preta_true,
              "overlap_preta_true" = overlap_preta_true
              )
         )
  },error=function(e){
    return(NA)
  }
  )
  
}

clusterExport(cl,"X1")

clusterExport(cl,"NT")

clusterExport(cl,"N.sim") 

clusterExport(cl, "models")

clusterExport(cl, "%!in%")

clusterExport(cl, "averageoverlap")

clusterExport(cl, "alternativenma")

clusterExport(cl, "nma.krahn.output")

clusterEvalQ(cl, {library(mtrank)})

clusterEvalQ(cl, {library(ircor)})

clusterEvalQ(cl, {library(tidyverse)})

clusterEvalQ(cl, {library(netmeta)})

clusterEvalQ(cl, {library(rlist)})

l <- parLapply(cl,1:N.sim, function(x) models(X1[[x]]))

l <- l[!is.na(l)]

N.sim <- length(l)

### Ability based metric vs true

kendall_PL_true=c()
for (i in 1:N.sim){kendall_PL_true=c(kendall_PL_true, l[[i]]$kendall_PL_true)}
mean(kendall_PL_true)

spearman_PL_true=c()
for (i in 1:N.sim){spearman_PL_true=c(spearman_PL_true, l[[i]]$spearman_PL_true)}
mean(spearman_PL_true)

overlap_PL_true=c()
for (i in 1:N.sim){overlap_PL_true=c(overlap_PL_true, l[[i]]$overlap_PL_true)}
mean(overlap_PL_true)

###### P-scores (random) vs true

kendall_prandom_true=c()
for (i in 1:N.sim){kendall_prandom_true=c(kendall_prandom_true, l[[i]]$kendall_prandom_true)}
mean(kendall_prandom_true)

spearman_prandom_true=c()
for (i in 1:N.sim){spearman_prandom_true=c(spearman_prandom_true, l[[i]]$spearman_prandom_true)}
mean(spearman_prandom_true)

overlap_prandom_true=c()
for (i in 1:N.sim){overlap_prandom_true=c(overlap_prandom_true, l[[i]]$overlap_prandom_true)}
mean(overlap_prandom_true)

###### P-scores (common) vs true

kendall_pcommon_true=c()
for (i in 1:N.sim){kendall_pcommon_true=c(kendall_pcommon_true, l[[i]]$kendall_pcommon_true)}
mean(kendall_pcommon_true)

spearman_pcommon_true=c()
for (i in 1:N.sim){spearman_pcommon_true=c(spearman_pcommon_true, l[[i]]$spearman_pcommon_true)}
mean(spearman_pcommon_true)

overlap_pcommon_true=c()
for (i in 1:N.sim){overlap_pcommon_true=c(overlap_pcommon_true, l[[i]]$overlap_pcommon_true)}
mean(overlap_pcommon_true)

###### pBV random vs true

kendall_pBV_random_true <- c()
for (i in 1:N.sim){kendall_pBV_random_true=c(kendall_pBV_random_true, l[[i]]$kendall_pBV_random_true)}
mean(kendall_pBV_random_true)

spearman_pBV_random_true  <- c()
for (i in 1:N.sim){spearman_pBV_random_true = c(spearman_pBV_random_true, l[[i]]$spearman_pBV_random_true)}
mean(spearman_pBV_random_true)

overlap_pBV_random_true=c()
for (i in 1:N.sim){overlap_pBV_random_true=c(overlap_pBV_random_true, l[[i]]$overlap_pBV_random_true)}
mean(overlap_pBV_random_true)

###### PRETA vs true 

kendall_preta_true=c()
for (i in 1:N.sim){kendall_preta_true=c(kendall_preta_true, l[[i]]$kendall_preta_true)}
mean(kendall_preta_true)

spearman_preta_true=c()
for (i in 1:N.sim){spearman_preta_true=c(spearman_preta_true, l[[i]]$spearman_preta_true)}
mean(spearman_preta_true)

overlap_preta_true=c()
for (i in 1:N.sim){overlap_preta_true=c(overlap_preta_true, l[[i]]$overlap_preta_true)}
mean(overlap_preta_true)

## Models vs true (kendall)

models_true_kendall <- cbind.data.frame(kendall_PL_true,
                                        kendall_prandom_true,
                                        kendall_pcommon_true,
                                        kendall_preta_true,
                                        kendall_pBV_random_true
                                        )

names(models_true_kendall) <- c("Ability based metric",
                                "P-score (random)",
                                "P-score (common)",
                                "PRETA (random)",
                                "pBV (random)"
                                )

avg_models_true_kendall <- cbind.data.frame(colMeans(models_true_kendall))

names(avg_models_true_kendall) <- c("avg_kendall") 

### print results for Kendall coefficient

avg_models_true_kendall

##########################################

## Models vs true (spearman)

models_true_spearman <- cbind.data.frame(spearman_PL_true,
                                         spearman_prandom_true,
                                         spearman_pcommon_true,
                                         spearman_preta_true,
                                         spearman_pBV_random_true
                                         )

names(models_true_spearman) <- c("Ability based metric",
                                "P-score (random)",
                                "P-score (common)",
                                "PRETA (random)",
                                "pBV (random)"
)



avg_models_true_spearman <- cbind.data.frame(colMeans(models_true_spearman))

names(avg_models_true_spearman) <- c("avg_spearman") 

### print results for Spearman coefficient

avg_models_true_spearman

##########################################


## Models vs true (overlap)

models_true_overlap <- cbind.data.frame(overlap_PL_true,
                                        overlap_prandom_true,
                                        overlap_pcommon_true,
                                        overlap_preta_true,
                                        overlap_pBV_random_true
                                        )

names(models_true_overlap) <- c("Ability based metric",
                                "P-score (random)",
                                "P-score (common)",
                                "PRETA (random)",
                                "pBV (random)"
)

avg_models_true_overlap <- cbind.data.frame(colMeans(models_true_overlap))

names(avg_models_true_overlap) <- c("avg_overlap") 

### print results for overlap

avg_models_true_overlap
