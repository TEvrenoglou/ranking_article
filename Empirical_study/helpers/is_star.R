is_star <- function(data){
  
chk1 <- length(unique(data$treat1))

chk2 <- length(unique(data$treat2))
  
res <- c()

if((chk1==1) || (chk2==1)){
  
res <- TRUE  
  
}else{
  res <- FALSE
}

return(res)
  
}

# count <- list()
# 
# for(i in 1:length(mod_RR)){
#   
#   
# count[[i]] <- is_star(mod_RR[[i]]$data)  
#   
# }
# 
# 
# stars <- table(unlist(count))
