pscore_graph_t<-function(x,beta_,type){
  
  pscores_applied=matrix(rep(0,length(beta_)*length(prep(x)$comm)),length(beta_),length(prep(x)$comm))
  # beta_<-beta_function(type,beta_)
  colnames(pscores_applied)<-prep(x)$comm
  for (i in 1:length(beta_)){
    pscores_applied[i,]<-pscores(prep(x)$outcomes,prep(x)$var.outcomes,1,-beta_[i],type,prep(x)$comm)
    par(oma=c(0, 0, 0, 0))
    plot(0,0,ylim=c(0,1),type="l",xlim=c(min(beta_),max(beta_)),ylab="p-score",xlab="CIV",main="")
    
    df<-data.frame(pscores_applied[1,],prep(x)$comm)
    names(df)<-c("pscore","interv")
    df_sort<-df[order(df$pscore),]
    best<-(tail(df_sort,n=5))$interv
    
    pscores_applied_1<-pscores_applied[,best]
    pscores_applied_2<-pscores_applied[,-which(colnames(pscores_applied) %in% best)]}
  
  for (j in 1:length(best)) {
    colr_b<-c("red","green","blue","orange","purple")
    lines(beta_,rowMeans(pscores_applied_1),type='b')
    lines(beta_,pscores_applied_1[1:length(beta_),j],col=colr_b[j]) 
  }
  legend("topright", legend=best,lwd = 2, col = colr_b,title="Best Interventions")
  
  for (j in length(best):length(prep(x)$comm)) {
    colr<-c(rep("black",length(prep(x)$comm)-5))
    lines(beta_,pscores_applied_2[1:length(beta_),j],col=colr[j]) 
  }
}