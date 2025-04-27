### Load the necessary libraries, all of which are available in CRAN. 
library(mtrank)
library(netmeta)
library(meta)
library(metafor)
library(tidyverse)
library(nmajags)
library(R2jags)


### source code to get ranking based on the PReTA approach
### Code obtained from the GitHub page: https://github.com/esm-ispm-unibe-ch/alternativenma/tree/master/R

source("./helpers/nma.krahn.output.R")

source("./helpers/alternative_nma.R")

source("./helpers/netmetaranks.R")


### source code to get ranking based on the adjusted according to MCID P-scores
### Code obtained from the GitHub page: https://github.com/DimitrisMavridis/RankingNMA/blob/master/extendedP-scores

source("./helpers/Pscores_function.R")
source("./helpers/league_table.R")
source("./helpers/league_table_var.R")
source("./helpers/intersect2.R")
source("./helpers/Prepare_Multi.R")
source("./helpers/Prepare_Single.R")
source("./helpers/Prepare_function.R")
source("./helpers/pscrs.R")
source("./helpers/pscore_graph.R")


### source code to get ranking based on the P(best) Bayesian approach
source("./helpers/Bayesian_NMA_model.R")
source("./helpers/get_results.R")


### get the 'diabetes' dataset from the mtrank package

data("diabetes")


### RANK TREATMENTS BASED ON P-scores ###

p <- pairwise(data = diabetes,
              event =  r,
              studlab = study,
              treat = t,
              n = n,
              sm = "OR"
)

## fit a NMA model
mod_netmeta <- netmeta(p,reference.group = "Placebo")

## visualize NMA estimates 
forest(mod_netmeta,
       overall.hetstat = TRUE, addrows = 0, calcwidth.hetstat = TRUE,
       print.I2 = FALSE, print.tau = TRUE, print.pval.Q = FALSE,
       digits.tau = 2,
       #
       drop.reference.group = TRUE,
       #
       label.left = "Favors vortioxetine",
       label.right = "Favors other treatments",
       #
       header.line = TRUE, spacing = 1.5,
       #
       file = "forest-diabetes-netmeta.pdf")

## get ranking based on P-scores
netrank(mod_netmeta,small.values = "desirable")


### RANK TREATMENTS BASED ON THE PROPOSED APPROACH ###

### Define a tcc using 1.20 as the MCID

# help(tcc)

ranks <- tcc(mod_netmeta,
             mcid = 1.20,
             small.values = "desirable"
)


## Fit the model and get the ability estimates

# help(mtrank)

mod_ability <- mtrank(ranks)

## Extract ability estimates
mod_ability$estimates

## Extract probability that each treatment is ranked first
mod_ability$probabilities

## Vizualize the ability estimates
forest(mod_ability, spacing = 1.5,
       file = "forest-diabetes-mtrank.pdf")

## Calculate pairwise probabilities

#help(paired_pref) 

fitted(mod_ability,treat1 = "ARB",treat2 = "Diuretic",type = "all")

fitted(mod_ability,treat1 = "ARB",treat2 = "ACE",type = "all")

### RANK TREATMENTS BASED ON THE PReTA approach ###

mod_preta <- alternativenma(mod_netmeta,small.values = "good")

## get ranking results
res_preta <- cbind.data.frame("treat"=row.names(mod_preta$averages),"PReTA_ranking"=mod_preta$averages$Pscoreaverage)

res_preta <- res_preta %>% 
  arrange(desc(PReTA_ranking))

res_preta

### RANK TREATMENTS BASED ON THE ADJUSTED ACCORDING TO THE MCID P-scores ###

net1 <- list(mod_netmeta)

pscores_mcid <- p_scores(net1,log(1.20),NULL,"H")

res_pscores_mcid <- cbind.data.frame("treat"=names(pscores_mcid),"P_score"=unname(pscores_mcid))

res_pscores_mcid  <- res_pscores_mcid  %>% 
  arrange(desc(P_score))

res_pscores_mcid

### RANK TREATMENTS BASED ON THE BAYESIAN P(best) approach ###

run.data = long2jags(
  studlab = study,
  treat = t,
  event = r,
  n = n,
  data = diabetes
)

drug_names <- levels(as.factor(diabetes$t))

run.list = list(
  ns = run.data$k,
  nt = run.data$n,
  na = run.data$n.k,
  t = run.data$T,
  r = run.data$E,
  n = run.data$N,
  ref = which(drug_names=="Placebo"),
  outcome=0 ## small.values are desirable
)

## Fit the Bayesian model

mod_bayesian = jags(
  data = run.list,
  inits = NULL,
  parameters.to.save = c(
    "logOR.ref",
    "tau",
    "effectiveness",
    "most.effective"
  ),
  n.chains = 2,
  n.iter = 50000,
  n.burnin = 10000,
  DIC = T,
  model.file = "NMA.model.txt"
)

## get results

results_Bayesian = as.data.frame(mod_bayesian$BUGSoutput$summary)

results_Bayesian$ind = as.character(rownames(results_Bayesian))

## extract P(best)

p_best = results_Bayesian %>%
  filter(grepl("most.effective", ind))

p_best <- data.frame("treatment"=drug_names,
                     "pBV"=p_best[,1])

p_best <- p_best %>% 
  arrange(desc(pBV))

p_best
