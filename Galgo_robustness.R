

load("/home/mguerrero/Genetic_alg/App_FINAL/Data/RNA_BRCA.rda")
downl=names(esets)

##Galgo Hyperparameters

population= 300                   #Number of individuals to evaluate
generations=500                 #Number of generations
nCV=5                             #Number of crossvalidations for function "crossvalidation"
GPU= TRUE                         # to use gpuR
TournamentSize=2
period=5475


set.seed(2003)
SET= sample(1:dim(esets[["metabric"]])[2],dim(esets[["metabric"]])[2]*0.75,replace=FALSE)

trainSet= esets[["metabric"]][,SET]
validationSet= esets[["metabric"]][,-SET]

prob_matrix= exprs(trainSet)
clinical=pData(trainSet)
OS=Surv(time=as.numeric(as.character(clinical$time)),event=clinical$status)
chrom_length= nrow(prob_matrix)   #length of chromosome

prob_matrixV= exprs(validationSet)
clinicalV= pData(validationSet)
OSv=Surv(time=as.numeric(as.character(clinicalV$time)),event=clinicalV$status)

#Check order
identical(colnames(prob_matrix),rownames(clinical))

#RUN alg for 30 runs

source("./Functions/functions.R")

for(h in 1:30){
  dir.create(paste0("./Results/Results_",h),recursive=TRUE)
  resultdir=paste0("./Results/Results_",h)
  source("./Functions/geneticalg.R")
}


#RUNS comparison



RESS= list()
for(j in 1:30){
  RES=data.frame(solution=as.character(),k=as.numeric(),ngenes=as.numeric(),trainSil=as.numeric(),testSil=as.numeric(),trainC=as.numeric(),testC=as.numeric(),ps=as.numeric(),stringsAsFactors = FALSE)
  
  resultdir=paste0("./Results/Results","_",j,"/")
  file= paste0("results",500,".rda")
  load(paste0(resultdir,file))
  X1=output$Solutions
  PARETO=output$ParetoFront
  dots=do.call(rbind,PARETO)
  
  rownames(X1)=paste("result",1:nrow(X1),sep=".")
  
  for(i in 1:sum(X1[,"rnkIndex"]==1)){
    
    RESULT=names(which(X1[order(X1[,(chrom_length+2)]),"rnkIndex"]==1))[i]
    R=as.logical(X1[RESULT,1:chrom_length])
    k=X1[RESULT,"k"]
    D=mydist(prob_matrix[R,])
    
    hsp_class=mycluster2(D,k)$cluster
    C=kcentroid(prob_matrix[R,],hsp_class)
    
    hsp_class= as.factor(classify(prob_matrix[R,],C,method=distancetype))
    hsp_classdf=data.frame(hsp_class=as.factor(hsp_class))
    
    mysurv=OS
    coxsimple=coxph(mysurv~hsp_class,data=hsp_classdf)
    a=concordance.index(predict(coxsimple),surv.time=OS[,1],surv.event=OS[,2],outx=FALSE)$c.index
    
    c=mean(silhouette(as.numeric(hsp_class),D)[,3])
    
    #Validation Set
    hsp_class= as.factor(classify(prob_matrixV[R,],C,method=distancetype))
    hsp_classdf=data.frame(hsp_class=as.factor(hsp_class))

    D2=mydist(prob_matrixV[R,])
    m=mean(silhouette(as.numeric(hsp_class),D2)[,3])
    
    hsp_class2=mycluster2(D2,k)$cluster
    
    ps=ps.cluster(hsp_class,hsp_class2)$ps
   
    d=concordance.index(predict(coxsimple,hsp_classdf),surv.time=OSv[,1],surv.event=OSv[,2],outx=FALSE)$c.index
    row= c(RESULT,k,sum(R),c,m,a,d,ps)
    RES[nrow(RES)+1,]=row
    
  }
  RESS[[j]]=RES
  
}

resAll=function(RESS){
  XX=lapply(RESS, function(x) { x["solution"] <- NULL; x })
  XX=lapply(XX, data.matrix)
  nsol=sapply(XX,nrow)
  all=sapply(XX,colMeans)
  ress= rbind(nsol,all)
  return(ress)
}

resBOX=function(RESS,feat){
  XX=lapply(RESS, function(x) { x["solution"] <- NULL; x })
  XX=lapply(XX, data.matrix)
  
  require(reshape2)
  df <- melt(XX)
  # plot using base R boxplot function
  boxplot(value ~ L1, data = df[df$Var2==feat,])
}
BOX=resAll(RESS)
resBOX(RESS,"trainC")

RESSsave=RESS
for(i in 1:length(RESS)){
  RESS[[i]]$run=i
  RESS[[i]]$weigth=1/nrow(RESS[[i]])
  
}
ALL=do.call(rbind,RESS)

mean(table(ALL[,"run"]))
min(table(ALL[,"run"]))
max(table(ALL[,"run"]))

min(ALL[,"trainC"])
max(ALL[,"trainC"])
weighted.mean(as.numeric(ALL[,"trainC"]),ALL[,"weigth"])

min(ALL[,"testC"])
max(ALL[,"testC"])
weighted.mean(as.numeric(ALL[,"testC"]),ALL[,"weigth"])

weighted.mean(as.numeric(ALL[,"k"]),ALL[,"weigth"])

table(ALL[,"run"],ALL[,"k"])

min(ALL[,"trainSil"])
max(ALL[,"trainSil"])
weighted.mean(as.numeric(ALL[,"trainSil"]),ALL[,"weigth"])

min(ALL[,"testSil"])
max(ALL[,"testSil"])
weighted.mean(as.numeric(ALL[,"testSil"]),ALL[,"weigth"])


cor(as.numeric(ALL[,"k"]),as.numeric(ALL[,"testSil"]),method="spearman)

cor(as.numeric(ALL[,"trainC"]),as.numeric(ALL[,"trainSil"]),method="spearman)

cor(as.numeric(ALL[,"testC"]),as.numeric(ALL[,"testSil"]),method="spearman)

cor(as.numeric(ALL[,"testC"]),as.numeric(ALL[,"trainC"]),method="spearman)
cor(as.numeric(ALL[,"testSil"]),as.numeric(ALL[,"trainSil"]),method="spearman)

library(ppcor)
pcor(data.frame(as.numeric(ALL$trainSil),as.numeric(ALL$trainC),as.numeric(ALL$k)),method="spearman")
cor(as.numeric(ALL$trainSil),as.numeric(ALL$trainC),method="spearman")

pcor(data.frame(as.numeric(ALL$testSil),as.numeric(ALL$testC),as.numeric(ALL$k)),method="spearman")
cor(as.numeric(ALL$testSil),as.numeric(ALL$testC),method="spearman")


pcor(data.frame(as.numeric(ALL$trainSil),as.numeric(ALL$ps),as.numeric(ALL$k)),method="spearman")

library(reshape2)
ALL2= melt(ALL[,c("trainSil","testSil","trainC","testC","ps","run")],id=c("run"))
ALL2$run=as.factor(as.numeric(ALL2$run))
ALL2$value=as.numeric(ALL2$value)

library(ggthemes)
ggplot(ALL2, aes(x=run, y=value, fill=variable)) + 
    geom_boxplot(outlier.size=0.5, outlier.shape=NA) +theme_hc() + scale_fill_brewer(palette="Paired") + stat_summary(fun.y=median, geom="smooth", aes(group=variable),lwd=1,color=rep(brewer.pal(5, "Paired"),each=30) )+ylim(0,1)

ALL$ngenes=as.numeric(ALL$ngenes)
ALL$k=as.factor(as.numeric(ALL$k))
ALL$run=as.factor(as.numeric(ALL$run))
 ggplot(ALL, aes(x=run, y=ngenes,fill=k)) + 
    geom_dotplot(dotsize = 0.5,binaxis="y",binwidth = 2, stackratio =.2,stackdir="center") + theme_hc() + scale_fill_brewer(palette="Paired") +ylim(0,150)
  
ggplot(ALL, aes(x=run, y=ngenes,fill=k)) + theme_hc() + geom_jitter(width = 0.2,size=1.5,alpha=0.5,shape=19,aes(color=k))+ scale_colour_brewer(palette="Set1") +ylim(0,150)
