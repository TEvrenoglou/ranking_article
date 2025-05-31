
all_rankings <- function(data,models,swd){

ranks <- list()

small.values <- list()

small.values.preta <- list()

small.values.pciv <- list()

ests_mtrank <- list()

mod <- list()

prbs <- list()

pscore <- list()

pbv <- list()

preta <- list()

pciv <- list()

res_all <- list()

nmadb_id <- list()

id <- list()

id_old <- list()

id_mod <- list()

is.star <- list()

varTE.random <- list()

avg_var <- list()

range <- list()

nS <- list()

nT <- list()

nS_nT <- list()

tau_un <- list()

all.ties.id <- list()

for(i in 1:length(models)){
  
  small.values[[i]] <- models[[i]]$small.values
  
  small.values.preta[[i]] <- ifelse(small.values[[i]]=="desirable",
                                    "good",
                                    "bad"
  )
  
  small.values.pciv[[i]] <- ifelse(small.values[[i]]=="desirable",
                                   "H",
                                   "B"
  )
  

  ranks[[i]] <- tcc(models[[i]],
                          small.values = small.values[[i]],
                          swd = swd
                          )

if(isTRUE(ranks[[i]]$all.ties)){
  
 all.ties.id[[i]] <- i 
}
  
  tryCatch({

mod[[i]] <- mtrank(ranks[[i]])

prbs[[i]] <- mod[[i]]$probabilities

prbs[[i]] <- prbs[[i]] %>% 
  arrange(treatment)

ests_mtrank[[i]] <- mod[[i]]$estimates  

ests_mtrank[[i]] <- ests_mtrank[[i]] %>% 
  arrange(treatment)

ests_mtrank[[i]]$prbs <- prbs[[i]]$probability 

pscore[[i]] <- netrank(models[[i]])

pscore[[i]] <- data.frame("treatment"=names(pscore[[i]]$ranking.random),
                          "P_score" = unname(pscore[[i]]$ranking.random)
)

pscore[[i]] <- pscore[[i]] %>% 
  arrange(treatment)
                      
pbv[[i]] <- rankogram(models[[i]])

pbv[[i]] <- data.frame("treatment"=names(pbv[[i]]$ranking.matrix.random[,1]),
                       "pbv" = unname(pbv[[i]]$ranking.matrix.random[,1])
                       )
  
preta[[i]] <- alternativenma(models[[i]],small.values = small.values.preta[[i]])

## get PReTA ranking results
preta[[i]] <- cbind.data.frame("treatment"=row.names(preta[[i]]$averages),
                               "PReTA"=preta[[i]]$averages$Pscoreaverage)

preta[[i]] <- preta[[i]]%>% 
  arrange(treatment)

  
## get P_score (CIV) results
  
  pciv[[i]] <- p_scores(list(models[[i]]),log(swd),NULL,small.values.pciv[[i]])
  
  pciv[[i]] <- cbind.data.frame("treatment"=names(pciv[[i]]),"P_CIV"=unname(pciv[[i]]))
  
  pciv[[i]]  <- pciv[[i]]  %>% 
    arrange(treatment)
    
 
res_all[[i]] <- cbind.data.frame("id"=i,
                                 ests_mtrank[[i]],
                                 "P_score"=pscore[[i]]$P_score,
                                 "pBV"=pbv[[i]]$pbv,
                                 "PReTA"=preta[[i]]$PReTA,
                                 "P_score_swd"=pciv[[i]]$P_CIV,
                                 "tau" = models[[i]]$tau
                                 )  

res_all[[i]]$id_old <- unique(data[[i]]$id_old)

res_all[[i]]$nmadb_id <- unique(data[[i]]$dat_id)
  
res_all[[i]]$id_mod <- unique(data[[i]]$id_mod)

res_all[[i]]$is.star <- unique(data[[i]]$is.star)   

res_all[[i]] <- res_all[[i]] %>% 
  dplyr::select(nmadb_id,id,id_old,id_mod,is.star,treatment,tau,log_ability,se,lower,upper,prbs,
                P_score,pBV,PReTA,P_score_swd)

### get metrics 

nmadb_id[[i]] <- unique(res_all[[i]]$nmadb_id)

id[[i]] <- unique(res_all[[i]]$id)

id_old[[i]] <- unique(res_all[[i]]$id_old)

id_mod[[i]] <- unique(res_all[[i]]$id_mod)

is.star[[i]] <- unique(res_all[[i]]$is.star)

varTE.random[[i]] <- models[[i]]$seTE.random^2

varTE.random[[i]] <- varTE.random[[i]][lower.tri(varTE.random[[i]])]

avg_var[[i]] <- mean(varTE.random[[i]])

range[[i]] <- ((max(varTE.random[[i]]))-(min(varTE.random[[i]])))/(max(varTE.random[[i]]))

nS[[i]] <- sum(models[[i]]$data$n1+models[[i]]$data$n2)

nT[[i]] <- length(models[[i]]$trts)

nS_nT[[i]] <- nS[[i]]/nT[[i]]

tau_un[[i]] <- models[[i]]$tau
  
  },error=function(e){
  return(NA)
}
)
  }

res_all <- bind_rows(res_all)

attr(res_all,"nmadb_id") <- unlist(nmadb_id)

attr(res_all,"id") <- unlist(id)

attr(res_all,"id_old") <- unlist(id_old)

attr(res_all,"id_mod") <- unlist(id_mod)

attr(res_all,"is.star") <- unlist(is.star)

attr(res_all,"avg_var") <- unlist(avg_var)

attr(res_all,"range") <- unlist(range)

attr(res_all,"nS_nT") <- unlist(nS_nT)

attr(res_all,"tau") <- unlist(tau_un)

attr(res_all,"all.ties.id") <-  unlist(all.ties.id)

return(res_all)

}

