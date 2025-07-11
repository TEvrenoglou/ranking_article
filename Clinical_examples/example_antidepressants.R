### Load the necessary libraries, all of which are available in CRAN. 
library(mtrank)
library(netmeta)
library(tidyverse)


### source code to get ranking based on the PReTA approach
### Code obtained from the GitHub page: https://github.com/esm-ispm-unibe-ch/alternativenma/tree/master/R

source("./Clinical_examples/helpers/nma.krahn.output.R")
source("./Clinical_examples/helpers/alternative_nma.R")
source("./Clinical_examples/helpers/netmetaranks.R")


### source code to get ranking based on the adjusted according to SWD P-scores
### Code obtained from the GitHub page: https://github.com/DimitrisMavridis/RankingNMA/blob/master/extendedP-scores

source("./Clinical_examples/helpers/Pscores_function.R")
source("./Clinical_examples/helpers/league_table.R")
source("./Clinical_examples/helpers/league_table_var.R")
source("./Clinical_examples/helpers/intersect2.R")
source("./Clinical_examples/helpers/Prepare_Multi.R")
source("./Clinical_examples/helpers/Prepare_Single.R")
source("./Clinical_examples/helpers/Prepare_function.R")
source("./Clinical_examples/helpers/pscrs.R")
source("./Clinical_examples/helpers/pscore_graph.R")

### get the 'antidepressants' dataset from the mtrank package
data("antidepressants")

### RANK TREATMENTS BASED ON P-scores ###

p <- pairwise(data = antidepressants,
              event =  responders,
              studlab = studyid,
              treat = drug_name,
              n = ntotal,
              sm = "OR")
# fit a NMA model
mod_netmeta <- netmeta(p,reference.group = "Trazodone")

# visualize NMA estimates (main manuscript, Figure 3A)
forest(mod_netmeta,
       overall.hetstat = TRUE, addrows = 0, calcwidth.hetstat = TRUE,
       print.I2 = FALSE, print.tau = TRUE, print.pval.Q = FALSE,
       digits.tau = 2,
       #
       drop.reference.group = TRUE,
       #
       label.left = "Favors trazodone",
       label.right = "Favors other treatments",
       #
       header.line = TRUE, spacing = 1.5,
       #
       file = "Main_manuscript_Figure3A.pdf")

# get ranking based on P-scores visualize NMA estimates (main manuscript, Table 1)
netrank(mod_netmeta,small.values = "undesirable")

### RANK TREATMENTS BASED ON THE PROPOSED APPROACH ###

# Define a tcc using 1.20 as the SWD
# help(tcc)
ranks <- tcc(mod_netmeta,
             swd = 1.20,
             small.values = "undesirable"
)

# Visualize the TCC in terms of the basic parameters ("Treatments vs Trazodone")
# help(forest.tcc)
forest(ranks,label.right="Favors other",label.left = "Favors Trazodone")

# Fit the model and get the ability estimates
# help(mtrank)
mod_ability <- mtrank(ranks)

# Extract ability estimates
mod_ability$estimates

# Extract probability that each treatment is ranked first
mod_ability$probabilities

# Vizualize the ability estimates (main manuscript, Figure 3B)
forest(mod_ability, spacing = 1.5,
       file = "Main_manuscript_Figure3B.pdf")

# create linegraph for sensitivity analysis (Main Manuscript, Figure 4)
sensitivity <- linegraph(mod_ability,
                         swd = seq(1.10,1.50,by=0.10),
                         swd.ref = 1.20,
                         k = 6
)

sensitivity

# same colors as in Main Manuscript, Figure 4
# install.packages("ggsci")
# library(ggsci)
# sensitivity+scale_color_bmj()+ylab("Probability (Normalized ability)")

### RANK TREATMENTS BASED ON THE PReTA approach ###

mod_preta <- alternativenma(mod_netmeta,small.values = "bad") 
# get ranking results (main manuscript, Table 2)
res_preta <- cbind.data.frame("treat"=row.names(mod_preta$averages),"PReTA_ranking"=mod_preta$averages$Pscoreaverage)
res_preta <- res_preta %>% 
  arrange(desc(PReTA_ranking))
res_preta

### RANK TREATMENTS BASED ON THE ADJUSTED ACCORDING TO THE SWD P-scores ###

net1 <- list(mod_netmeta)
pscores_SWD <- p_scores(net1,log(1.20),NULL,"B")
res_pscores_SWD <- cbind.data.frame("treat"=names(pscores_SWD),"P_score"=unname(pscores_SWD))
res_pscores_SWD  <- res_pscores_SWD  %>% 
  arrange(desc(P_score))
# (main manuscript, Table 2)
res_pscores_SWD

### RANK TREATMENTS BASED ON THE FREQUENTIST pBV approach ###

r <- rankogram(mod_netmeta,small.values = "undesirable",nsim = 50000)
pBV <- as.data.frame(r$ranking.matrix.random[,1])
names(pBV) <- c("pBV")
pBV <- pBV %>% 
  arrange(desc(pBV)) %>% 
  mutate(pBV=round(pBV,digits = 2))
# (main manuscript, Table 2)
pBV

