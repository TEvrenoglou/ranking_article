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


### source code to get ranking based on the adjusted according to SWD P-scores
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

# fit a NMA model
mod_netmeta <- netmeta(p,reference.group = "Placebo")

# visualize NMA estimates 
forest(mod_netmeta,
       overall.hetstat = TRUE, addrows = 0, calcwidth.hetstat = TRUE,
       print.I2 = FALSE, print.tau = TRUE, print.pval.Q = FALSE,
       digits.tau = 2,
       #
       drop.reference.group = TRUE,
       #
       label.left = "Favors other treatments",
       label.right = "Favors Placebo",
       #
       header.line = TRUE, spacing = 1.5,
       #
       file = "forest-diabetes-netmeta.pdf")

# get ranking based on P-scores
netrank(mod_netmeta,small.values = "desirable")


### RANK TREATMENTS BASED ON THE PROPOSED APPROACH ###

### Define a tcc using 1.20 as the SWD

# help(tcc)

ranks <- tcc(mod_netmeta,
             mcid = 1.20,
             small.values = "desirable")

# Visualize the TCC in terms of the basic parameters ("Treatments vs Placebo")
# help(forest.tcc)
forest(ranks,label.right="Favors Placebo",label.left = "Favors other treatments")

# Fit the model and get the ability estimates
# help(mtrank)
mod_ability <- mtrank(ranks)

# Extract ability estimates
mod_ability$estimates

# Extract probability that each treatment is ranked first
mod_ability$probabilities

# Vizualize the ability estimates
forest(mod_ability, spacing = 1.5,
       file = "forest-diabetes-mtrank.pdf")

### RANK TREATMENTS BASED ON THE PReTA approach ###

mod_preta <- alternativenma(mod_netmeta,small.values = "good")
# get ranking results
res_preta <- cbind.data.frame("treat"=row.names(mod_preta$averages),"PReTA_ranking"=mod_preta$averages$Pscoreaverage)
res_preta <- res_preta %>% 
  arrange(desc(PReTA_ranking))
res_preta

### RANK TREATMENTS BASED ON THE ADJUSTED ACCORDING TO THE SWD P-scores ###

net1 <- list(mod_netmeta)
pscores_SWD <- p_scores(net1,log(1.20),NULL,"H")
res_pscores_SWD <- cbind.data.frame("treat"=names(pscores_SWD),"P_score"=unname(pscores_SWD))
res_pscores_SWD  <- res_pscores_SWD  %>% 
  arrange(desc(P_score))
res_pscores_SWD

### RANK TREATMENTS BASED ON THE FREQUENTIST pBV approach ###

r <- rankogram(mod_netmeta,small.values = "desirable",nsim = 50000)
pBV <- as.data.frame(r$ranking.matrix.random[,1])
names(pBV) <- c("pBV")
pBV <- pBV %>% 
  arrange(desc(pBV)) %>% 
  mutate(pBV=round(pBV,digits = 2))
pBV
