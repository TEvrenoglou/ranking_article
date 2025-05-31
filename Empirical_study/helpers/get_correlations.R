
get_correlations <- function(data){

ids <- unique(data$id)

nmadb_id <- attributes(data)$nmadb_id

id_old <- attributes(data)$id_old

id_mod <- attributes(data)$id_mod

is.star <- attributes(data)$is.star

avg_var <- attributes(data)$avg_var

range <- attributes(data)$range

nS_nT <- attributes(data)$nS_nT

tau <- attributes(data)$tau


dat1 <- list()

cor_ability_pscore <- list()

cor_ability_pbv <- list()

cor_ability_PReTA <- list()

cor_ability_pswd <- list()


cor_prbs_pscore <- list()

cor_prbs_pbv <- list()

cor_prbs_PReTA <- list()

cor_prbs_pswd <- list()


cor_pscore_pbv <- list() 

cor_pscore_pswd <- list()

cor_pscore_PReTA <- list()


cor_pbv_PReTA <- list()

cor_pbv_pswd <- list()

cor_PReTA_pswd <- list()


for(i in 1:length(ids)){

dat1[[i]] <- data %>%
  filter(id==ids[i])
  
cor_ability_pscore[[i]] <- cor(dat1[[i]]$log_ability,dat1[[i]]$P_score)

cor_ability_pbv[[i]] <- cor(dat1[[i]]$log_ability,dat1[[i]]$pBV)

cor_ability_PReTA[[i]] <- cor(dat1[[i]]$log_ability,dat1[[i]]$PReTA)

cor_ability_pswd[[i]] <- cor(dat1[[i]]$log_ability,dat1[[i]]$P_score_swd)


cor_prbs_pbv[[i]] <- cor(dat1[[i]]$prbs,dat1[[i]]$pBV)

cor_prbs_pscore[[i]] <- cor(dat1[[i]]$prbs,dat1[[i]]$P_score)

cor_prbs_PReTA[[i]] <- cor(dat1[[i]]$prbs,dat1[[i]]$PReTA)

cor_prbs_pswd[[i]] <- cor(dat1[[i]]$prbs,dat1[[i]]$P_score_swd)



cor_pscore_pbv[[i]] <- cor(dat1[[i]]$P_score,dat1[[i]]$pBV)

cor_pscore_PReTA[[i]]<- cor(dat1[[i]]$P_score,dat1[[i]]$PReTA)

cor_pscore_pswd[[i]]<- cor(dat1[[i]]$P_score,dat1[[i]]$P_score_swd)


cor_pbv_PReTA[[i]]<- cor(dat1[[i]]$pBV,dat1[[i]]$PReTA)

cor_pbv_pswd[[i]]<- cor(dat1[[i]]$pBV,dat1[[i]]$P_score_swd)

cor_PReTA_pswd[[i]] <- cor(dat1[[i]]$PReTA,dat1[[i]]$P_score_swd)

}

cor_ability_pscore <- unlist(cor_ability_pscore)

cor_ability_pbv <- unlist(cor_ability_pbv)

cor_ability_PReTA <- unlist(cor_ability_PReTA)

cor_ability_pswd <- unlist(cor_ability_pswd)

cor_prbs_pbv <- unlist(cor_prbs_pbv)

cor_prbs_pscore <- unlist(cor_prbs_pscore)

cor_prbs_PReTA <- unlist(cor_prbs_PReTA)

cor_prbs_pswd <- unlist(cor_prbs_pswd)


cor_pscore_pbv <- unlist(cor_pscore_pbv)

cor_pscore_PReTA <- unlist(cor_pscore_PReTA)

cor_pscore_pswd <- unlist(cor_pscore_pswd)


cor_pbv_PReTA <- unlist(cor_pbv_PReTA)

cor_pbv_pswd <- unlist(cor_pbv_pswd)

cor_PReTA_pswd <- unlist(cor_PReTA_pswd)

all_cors <- cbind.data.frame("nmadb_id"=nmadb_id,"id"=ids,"id_old"=id_old,"id_mod"=id_mod,"is.star"=is.star,
                             cor_ability_pscore,cor_ability_pbv,cor_ability_PReTA,cor_ability_pswd,
                             cor_prbs_pscore,cor_prbs_pbv,cor_prbs_PReTA,cor_prbs_pswd,
                             cor_pscore_pbv,cor_pscore_PReTA,cor_pscore_pswd,
                             cor_pbv_PReTA,cor_pbv_pswd,
                             cor_PReTA_pswd,
                             "tau"=tau,
                             "avg_var"=avg_var,
                             "range"=range,
                             "nS_nT"=nS_nT
                             )

return(all_cors)

}

