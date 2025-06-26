library(tidyverse)
library(mtrank)
### Functions to get PReTA and adjusted according to P-score
source("./Empirical_study/helpers/alternative_nma.R")
source("./Empirical_study/helpers/nma.krahn.output.R")
source("./Empirical_study/helpers/Pscores_function.R")
source("./Empirical_study/helpers/league_table.R")
source("./Empirical_study/helpers/league_table_var.R")
source("./Empirical_study/helpers/intersect2.R")
source("./Empirical_study/helpers/Prepare_Multi.R")
source("./Empirical_study/helpers/Prepare_Single.R")
source("./Empirical_study/helpers/Prepare_function.R")
source("./Empirical_study/helpers/pscrs.R")
source("./Empirical_study/helpers/pscore_graph.R")
source("./Empirical_study/helpers/is_star.R")
# a list with networks from nmadb
nmadb <- readRDS("./Empirical_study/NMAs_nmadb/nmadb.rds")
# a list with NMAs in terms of Risk Ratio obtained from nmadb's networks
mod_RR <- readRDS("./Empirical_study/NMAs_nmadb/mod_RR.rds")
# a function to calculate rankings across methods 
source("./Empirical_study/helpers/all_rankings.R")
# a function to calculate correlation coeffiecients across methods 
source("./Empirical_study/helpers/get_correlations.R")

# get all ranks (21 warnings are expected from the 21 networks where only ties were identified)
# runs in approximately 4 minutes
all_ranks <- all_rankings(nmadb,mod_RR,1.25)
# get correlations across methods
cors <- get_correlations(all_ranks)

# Save the Results in the folder "Results". 
# You may skip this as Results are already written in the folder
res_save <- paste(getwd(),"/Empirical_study/Results",sep = "")
write.csv(all_ranks,paste(res_save,"/all_ranks.csv",sep = ""),row.names = FALSE)
write.csv(cors,paste(res_save,"/all_cors.csv",sep = ""),row.names = FALSE)

# Calculate correlations as in Table 3 (main manuscript)
# Summary results are presented in terms of median correlation and interquartile range

# Ability based metric vs rest methods
ab_pscr <- round(quantile(cors$cor_ability_pscore,probs = c(0.25,0.5,0.75)),digits = 3)
ab_pbv <- round(quantile(cors$cor_ability_pbv,probs = c(0.25,0.5,0.75)),digits = 3)
ab_preta <- round(quantile(cors$cor_ability_PReTA,probs = c(0.25,0.5,0.75)),digits = 3)
ab_pswd <- round(quantile(cors$cor_ability_pswd,probs = c(0.25,0.5,0.75)),digits = 3)

# P-score vs rest methods
pscr_pbv <- round(quantile(cors$cor_pscore_pbv,probs = c(0.25,0.5,0.75)),digits = 3)
pscr_preta <- round(quantile(cors$cor_pscore_PReTA,probs = c(0.25,0.5,0.75)),digits = 3)
pscr_pswd <- round(quantile(cors$cor_pscore_pswd,probs = c(0.25,0.5,0.75)),digits = 3)

# pBV vs rest methods
pbv_preta <- round(quantile(cors$cor_pbv_PReTA,probs = c(0.25,0.5,0.75)),digits = 3)
pbv_pswd <- round(quantile(cors$cor_pbv_pswd,probs = c(0.25,0.5,0.75)),digits = 3)

# PReTA vs rest methods
preta_pswd <- round(quantile(cors$cor_PReTA_pswd,probs = c(0.25,0.5,0.75)),digits = 3)

# Table 5 (main manuscript) small deviations due to re-sampling at the pbv method
res_cors <- rbind.data.frame(ab_pscr,ab_pbv,ab_preta,ab_pswd,
                             pscr_pbv,pscr_preta,pscr_pswd,
                             pbv_preta,pbv_pswd,
                             preta_pswd)

names(res_cors) <- c("25%","median","75%")

res_cors <- res_cors %>% 
  dplyr::select(median,'25%','75%')

row.names(res_cors) <- c("Ability_based_metric_vs_P-score",
                     "Ability_based_metric_vs_pBV",
                     "Ability_based_metric_vs_PReTA",
                     "Ability_based_metric_vs_P-score_(SWD)",
                     "P-score_vs_pBv",
                     "P-score_vs_PReTA",
                     "P-score_vs_P-score_(SWD)",
                     "pBV_vs_PReTA",
                     "pBV_vs_P-score_(SWD)",
                     "PReTA_vs_P-score_(SWD)")
                     
res_cors

write.csv(res_cors,paste(res_save,"/Table3_main_manuscript.csv",sep = ""),row.names = TRUE)

# Reproduce Figure 7 of the main manuscript

filename_var <- "./Empirical_study/Results/Figure_7A_main_manuscript.tif"
tiff(filename = filename_var,res = 300,width = 2500,height = 2500)
par(mfrow = c(2, 2), oma = c(4, 5, 2, 1))  # Set outer margins: bottom, left, top, right
plot(log(cors$avg_var),cors$cor_ability_pscore,xlab = "",ylab = "",main = "Ability based metric vs P-scores")
lines(smooth.spline(log(cors$avg_var), cors$cor_ability_pscore,df=5),lwd=2,col="purple")
plot(log(cors$avg_var),cors$cor_ability_PReTA,xlab = "",ylab = "",main = "Ability based metric vs PReTA")
lines(smooth.spline(log(cors$avg_var), cors$cor_ability_PReTA,df=5),lwd=2,col="purple")
plot(log(cors$avg_var),cors$cor_ability_pbv,xlab = "",ylab = "",main = "Ability based metric vs pBV")
lines(smooth.spline(log(cors$avg_var), cors$cor_ability_pbv,df=5),lwd=2,col="purple")
plot(log(cors$avg_var),cors$cor_ability_pswd,xlab = "",ylab = "",main = "Ability based metric vs P-scores (SWD)")
lines(smooth.spline(log(cors$avg_var), cors$cor_ability_pswd,df=5),lwd=2,col="purple")
mtext("log(average variance)", side = 1, outer = TRUE, line = 1, cex = 1.5)
mtext("Pearson's correlation coefficient", side = 2, outer = TRUE, line = 1, cex = 1.5)
dev.off()

filename_range <- "./Empirical_study/Results/Figure_7B_main_manuscript.tif"
tiff(filename = filename_range,res = 300,width = 2500,height = 2500)
par(mfrow = c(2, 2), oma = c(4, 5, 2, 1))  # Set outer margins: bottom, left, top, right
plot(log(log(1/cors$range)),cors$cor_ability_pscore,xlab = "",ylab = "",main = "Ability based metric vs P-scores")
lines(smooth.spline(log(log(1/cors$range)), cors$cor_ability_pscore,df=5),lwd=2,col="purple")
plot(log(log(1/cors$range)),cors$cor_ability_PReTA,xlab = "",ylab = "",main = "Ability based metric vs PReTA")
lines(smooth.spline(log(log(1/cors$range)), cors$cor_ability_PReTA,df=5),lwd=2,col="purple")
plot(log(log(1/cors$range)),cors$cor_ability_pbv,xlab = "",ylab = "",main = "Ability based metric vs pBV")
lines(smooth.spline(log(log(1/cors$range)), cors$cor_ability_pbv,df=5),lwd=2,col="purple")
plot(log(log(1/cors$range)),cors$cor_ability_pswd,xlab = "",ylab = "",main = "Ability based metric vs P-scores (SWD)")
lines(smooth.spline(log(log(1/cors$range)), cors$cor_ability_pswd,df=5),lwd=2,col="purple")
mtext("log(log(1/relative range variance))", side = 1, outer = TRUE, line = 1, cex = 1.5)
mtext("Pearson's correlation coefficient", side = 2, outer = TRUE, line = 1, cex = 1.5)
dev.off()

# Analysis of 21 networks with only ties
# a list of 21 NMAs containing only ties when the TCC is applied to them
ties <- readRDS("./Empirical_study/NMAs_nmadb/modID.ties.rds")

mod.ties <- small.vals <- r.ties <- TE <- seTE <- TE_des <- seTE_des <- varTE_des <-
  TE_undes <- seTE_undes <- varTE_undes <- right <- left <- star <- vector("list")

for (i in 1:length(ties)){
  
  mod.ties[[i]] <- mod_RR[[ties[i]]]
  
  star[[i]] <- is_star(mod.ties[[i]]$data)
  
  small.vals[[i]] <- mod.ties[[i]]$small.values
  
  r.ties[[i]] <- tcc(mod.ties[[i]],swd = 1.25,
                     small.values = small.vals[[i]])
  
  TE[[i]] <- mod.ties[[i]]$TE.random[lower.tri(mod.ties[[i]]$TE.random)]
  
  seTE[[i]] <- mod.ties[[i]]$seTE.random[lower.tri(mod.ties[[i]]$seTE.random)]
  
  if(small.vals[[i]] == "desirable"){
    
    TE_des[[i]] <- mod.ties[[i]]$TE.random[lower.tri(mod.ties[[i]]$TE.random)]
    
    seTE_des[[i]] <-  mod.ties[[i]]$seTE.random[lower.tri(mod.ties[[i]]$seTE.random)]
    
    varTE_des[[i]] <- seTE_des[[i]]^2
    
  }else{
    
    TE_des[[i]] <- NA
    
    seTE_des[[i]] <- NA
    
    varTE_des[[i]] <- NA
  }
  
  if(small.vals[[i]] == "undesirable"){
    
    TE_undes[[i]] <- mod.ties[[i]]$TE.random[lower.tri(mod.ties[[i]]$TE.random)]
    
    seTE_undes[[i]] <-  mod.ties[[i]]$seTE.random[lower.tri(mod.ties[[i]]$seTE.random)]
    
    varTE_des[[i]] <-  seTE_undes[[i]]^2
    
  }else{
    
    TE_undes[[i]] <- NA
    
    seTE_undes[[i]] <- NA
    
    varTE_undes[[i]] <- NA
  }
  
}

# treatment effects (log-scale) and standard errors in networks with harmful outcomes
TE_des <- Filter(function(x) !is.atomic(x) || length(x) != 1 || !is.na(x), TE_des)
seTE_des <- Filter(function(x) !is.atomic(x) || length(x) != 1 || !is.na(x), seTE_des)
# treatment effects (natural scale) in networks with harmful outcomes
exp_TE_des <- sapply(TE_des,exp)

# treatment effects (log-scale) and standard errors in networks with beneficial outcomes
TE_undes <- Filter(function(x) !is.atomic(x) || length(x) != 1 || !is.na(x), TE_undes)
seTE_undes <- Filter(function(x) !is.atomic(x) || length(x) != 1 || !is.na(x), seTE_undes)
# treatment effects (natural scale) in networks with beneficial outcomes
exp_TE_undes <- sapply(TE_undes,exp)

# combine all treatment effects (natural scale) and standard errors

TE_all <- c(unlist(exp_TE_des),unlist(exp_TE_undes))
se_all <- c(unlist(seTE_des),unlist(seTE_undes))

# Get summary results
# summary treatment effects harmful outcomes
summary(unlist(exp_TE_des))
# summary treatment effects beneficial outcomes
summary(unlist(exp_TE_undes))
# summary TEs for across all outcomes
summary(TE_all)
## check star networks
star <- unlist(star)
## TRUE = Star network,  FALSE = Not a star network
table(star) 

## Save forest plots in pdfs (Supplementary Figures 3-23, Panel (a))
filename <- "./Empirical_study/Results/Supplementary_Figures_3_23_panel_a.pdf"
pdf(file = filename,
    width = 7,
    height = 7)
for(i in 1:length(mod.ties)){
  right[[i]] <- ifelse(small.vals[[i]]=="undesirable","Favors other","Favors 1")
  left[[i]] <- ifelse(small.vals[[i]]=="undesirable","Favors 1","Favors other")
if(i %in% c(5,6,20)){
  forest(mod.ties[[i]],sortvar = "-Pscore",
  rightcols = c("effect","ci","Pscore"),
  label.left = left[[i]],
  label.right = right[[i]],
  drop.reference.group = FALSE,
  cid.below.null = 0.80,
  cid.above.null = 1.25,
  col.cid = "lightblue",
  header.line = "both",
  fill.equi = "lightblue",
  col.square = "black",
  col.square.lines = "black",
  squaresize = 0.7,
  xlim = c(0.7,1.35)
  )
  }else{
    forest(mod.ties[[i]],sortvar = "-Pscore",
           rightcols = c("effect","ci","Pscore"),
           label.left = left[[i]],
           label.right = right[[i]],
           drop.reference.group = FALSE,
           cid.below.null = 0.80,
           cid.above.null = 1.25,
           col.cid = "lightblue",
           header.line = "both",
           fill.equi = "lightblue",
           col.square = "black",
           col.square.lines = "black",
           squaresize = 0.7
    )
  }
}
dev.off()

## Save network graphs in pdfs (Supplementary Figures 3-23, Panel (b))
filename_net <- "./Empirical_study/Results/Supplementary_Figures_3_23_panel_b.pdf"
pdf(file = filename_net,
    width = 7,
    height = 7)

for(i in 1:length(mod.ties)){netgraph(mod.ties[[i]])}
dev.off()
