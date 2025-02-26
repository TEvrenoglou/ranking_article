dlong <- function(data){
  # prepare the data set for the long-format
  dat1 <- data %>%
    dplyr::select(studlab, t1, n1, r1) %>%
    rename("study" = "studlab", 
           "treat" = "t1", 
           "events" = "r1", 
           "n" = "n1")
  
  dat2 <- data %>%
    dplyr::select(studlab, t2, n2, r2) %>%
    rename("study" = "studlab", 
           "treat" = "t2", 
           "events" = "r2", 
           "n" = "n2")
  
  data <- rbind(dat1, dat2)
  data <- data %>%
    mutate(ID = paste(study, treat, sep ="")) %>%
    filter(!duplicated(ID)) %>%
    dplyr::select(-ID)
  
  data=data[order(data$study),]
  return(data)
}
