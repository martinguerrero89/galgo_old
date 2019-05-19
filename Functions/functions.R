##################
####Functions#####
##################


#Tidy METABRIC and OSLO data

expressiontidy2= function(Expression,ANNOT, CLINIC, SURV,BATCH){
  require(illuminaio)
  require(AnnotationDbi)
  require(limma)
  
  CLINIC$ER.Expr= gsub("\\+", "positive", CLINIC$ER.Expr)
  CLINIC$ER.Expr= gsub("\\-", "negative", CLINIC$ER.Expr)
  CLINIC$PR.Expr= gsub("\\+", "positive", CLINIC$PR.Expr)
  CLINIC$PR.Expr= gsub("\\-", "negative", CLINIC$PR.Expr)
  CLINIC$Her2.Expr= gsub("\\+", "positive", CLINIC$Her2.Expr)
  CLINIC$Her2.Expr= gsub("\\-", "negative", CLINIC$Her2.Expr)
  SURV[,2]= gsub("1", "deceased", SURV[,2])
  SURV[,2]= gsub("0", "living", SURV[,2])
  pData= data.frame(
  row.names= make.names(rownames(CLINIC)),
  sample_name= make.names(rownames(CLINIC)),
  alt_sample_name= NA,
  unique_patient_ID= NA,
  sample_type= "tumor",
  er= CLINIC$ER.Expr,
  pgr= CLINIC$PR.Expr,
  her2= CLINIC$Her2.Expr,
  tumor_size= CLINIC$size/10,
  T=NA,
  N=NA,
  age_at_initial_pathologic_diagnosis= CLINIC$age_at_diagnosis,
  grade=CLINIC$grade,
  dmfs_days=NA,
  dmfs_status=NA,
  days_to_tumor_recurrence=NA,
  recurrence_status=NA,
  days_to_death=SURV[,1],
  vital_status=SURV[,2],
  tissue=NA,
  treatment=CLINIC$Treatment,
  percent_normal_cells=NA,
  percent_stromal_cells=NA,
  percent_tumor_cells=NA,
  batch=BATCH,
  uncurated_author_metadata=NA,
  duplicates=NA)
  probes=rownames(Expression)
  ANNOT=ANNOT$probe[ANNOT$probe[,"Probe_Id"] %in% probes,]
  ANNOT= ANNOT[!is.na(ANNOT$Entrez_Gene_ID),]
  Expression=Expression[ANNOT[,"Probe_Id"],]
  fData= data.frame(probeset=ANNOT$Probe_Id,gene=ANNOT$Symbol, EntrezGene.ID=ANNOT$Entrez_Gene_ID,best_probe=TRUE,row.names=ANNOT$Probe_Id)
  colnames(Expression)=make.names(colnames(Expression))
  phenoData <- new("AnnotatedDataFrame",data=pData)
  featureData <- new("AnnotatedDataFrame",data=fData)
  ExpSet= ExpressionSet(assayData=Expression,phenoData=phenoData,featureData=featureData)
  return(ExpSet)
}


#Expand probesets for multiple genes for a single probe

expandProbesets <- function (eset, sep = "///")
{
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  y<- lapply(as.character(fData(eset)$gene), function(x) strsplit(x, sep))
  eset <- eset[order(sapply(x, length)), ]
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  y<- lapply(as.character(fData(eset)$gene), function(x) strsplit(x, sep))
  idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
  idy <- unlist(sapply(1:length(y), function(i) rep(i, length(y[[i]]))))
  xx <- !duplicated(unlist(x))
  idx <- idx[xx]
  idy <- idy[xx]
  x <- unlist(x)[xx]
  y <- unlist(y)[xx]
  
  eset <- eset[idx, ]
  featureNames(eset) <- x
  fData(eset)$EntrezGene.ID <- x
  fData(eset)$gene <- y
  return(eset)
  
}


#Distance type used in the algorithm
distancetype="pearson"          #Options are: "pearson","uncentered" (cosine),"spearman","euclidean"

distancetype=match.arg(distancetype,c("pearson","uncentered","euclidean","spearman"))
#Centered pearson distance
                       
if(distancetype=="pearson"){
  gpuDist=function(x){
    x=vclMatrix(x)
    mx=gpuR::colMeans(x)
    x2=gpuR::colMeans(x^2)
    mx2=mx^2
    sx=sqrt(x2-mx2)
    pDist=((gpuR::crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
    #pDist=gpuR::cov(x)/sx%o%sx
    return(as.dist(1-pDist))
  } 
}

#Spearman distance
if(distancetype=="spearman"){
  gpuDist=function(x){
    x=apply(x,2,rank)
    x=vclMatrix(x)
    mx=gpuR::colMeans(x)
    x2=gpuR::colMeans(x^2)
    mx2=mx^2
    sx=sqrt(x2-mx2)
    pDist=((gpuR::crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
    #pDist=gpuR::cov(x)/sx%o%sx
    return(as.dist(1-pDist))
  } 
}

if(distancetype=="uncentered"){ 
  gpuDist=function(m){
    mgpu=gpuR::vclMatrix(t(m))
    d2= gpuR::vclMatrix(1,ncol=1,nrow=dim(mgpu)[2])
    a1=gpuR::tcrossprod(mgpu,mgpu)
    a2= mgpu^2 %*% d2
    pDist= a1/ sqrt(gpuR::tcrossprod(a2,a2))
    return(as.dist(1-pDist))
  } 
  
}
if(distancetype=="euclidean"){
  gpuDist= function(m){
    mgpu=gpuR::vclMatrix(t(m))
    d=suppressWarnings(gpuR::distance(mgpu,mgpu,method="euclidean"))
    return(as.dist(d))
  }
}

#Distance functions

mydist=function(c) {gpuDist(c)}                     #c is expression matrix
mydist2=function(x) mydist(t(x))

#To use hierchical clustering 
myclust=function(d) {hclust(d,method="average")}    #d is a distance object

mycluster1 <- function(d, k) {                      #d is a distance object and k is number of clusters
  hc=cutree(myclust(d),k=k)
  return(list(cluster=hc))
}

#To use partition around medioids (PAM), This is the default for the the current aplication
mycluster2<- function(c,k){
  return(list(cluster=pam(c,k,cluster.only=TRUE,diss=TRUE,do.swap=TRUE,keep.diss=FALSE,keep.data=FALSE,pamonce=2)))
}


#cosine similarity
cosine= function(a,b){a %*% b / sqrt(a%*%a * b%*%b)}

#Function to calculate the centroids of different groups (classes)
kcentroid=function(data,class){
  L=list()
  c=unique(unlist(class))
  for(i in c){
    if(sum(unlist(class)==i)>1){
      x=rowMeans(data[,unlist(class)==i])
      L[[i]]=x
    }else{
      L[[i]]=data[,unlist(class)==i]
    }}
  L=t(do.call(rbind,L))
  return(L)
}

#Function to calculate the distance to the centroids

gpuCor=function(x,y){
  x=vclMatrix(as.matrix(x))
  mx=gpuR::colMeans(x)
  x2=gpuR::colMeans(x^2)
  mx2=mx^2
  sx=sqrt(x2-mx2)
  
  y=vclMatrix(as.matrix(y))
  my=gpuR::colMeans(y)
  y2=gpuR::colMeans(y^2)
  my2=my^2
  sy=sqrt(y2-my2)
  
  pDist=((gpuR::crossprod(x,y)/nrow(x))-(mx %o% my))/sx%o%sy
  return(pDist)
} 


corrF2= function(data,centroid){
  RankD=apply(data,2,rank)
  RankC=apply(centroid,2,rank)
  scores=gpuCor(RankC,RankD)
  return(as.matrix(scores))
}

#Given the distance to the centroids classify the samples
classify2= function(data,centroid,method=distancetype){
  if(method=="pearson"){
    R=cor(data,centroid,method="pearson")
    scores<-apply(R,1,which.max)
    return(scores)}
  if(method=="spearman"){
    R=cor(data,centroid,method="spearman")
    scores<-apply(R,1,which.max)
    return(scores)}
}

#gpuRversion
classify= function(data,centroid,method=distancetype){
  if(method=="pearson"){
    R=gpuCor(data,centroid)
    scores<-apply(R,1,which.max)
    return(scores)}
  
  if(method=="spearman"){
    R=corrF2(data,centroid)
    scores<-apply(R,2,which.max)
    return(scores)}
}


#harmonic mean tends to be robust with high outlayers (robust for overfitting)
#with high penalty on small values
hmean=function(a){1/mean(1/a)}


#Fitness by RMST (Restricted Mean Survival Time, https://www.bmj.com/content/357/bmj.j2250)
cDist= function(x){ #ad-hoc function, x is the RMST
  d=x[order(x)]     
  l=length(d)-1
  c=c(0,1)
  dif=as.numeric()
  for(i in 1:l){
    dif=c(dif,diff(d[c+i]))
  }
  return(hmean(dif)*l)
}

fitness=function(OS,clustclass){
  score=tryCatch({ 
    t=survival:::survmean(survfit(OS~clustclass),rmean=period)[[1]][,"*rmean"] #This function calculates the RMST (comes from package Survival) 
    cDist(t)},error=function(e)return(0)) #If cDist cannot be calculated, the difference is set to 0 (no difference between curves)
  return(score)
}

##Functions to vectorize crossvalidation

ArrayTrain=function(flds,Data){
  trainData=Data[,-flds]
}

ArrayTest=function(flds,Data){
  testData=Data[,flds]
}


subDist=function(flds,D){
  sub=subset(D,-flds)
}

alloc2=function(C){
  Ct=t(C)
  ord=galeShapley.marriageMarket(C,Ct)
  return(ord$engagements)
}


reord=function(C,ord){
  C=C[,ord]
}



#crossvalidation function based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3105299/
#Data is the expression matrix
#flds is a list with the indexes to partition the data in n different folds
#indv is the solution to test, namely, a binary vector to subset the genes to test
#The function returns two fitness: fit1 (the mean silhouette) and fit2 (the ad-hoc function to estimate the differences between curves)

crossvalidation= function(Data, flds,indv,k,OS){
  Data= Data[indv,] 
  D=mydist(Data)
  TrainA= sapply(flds,ArrayTrain,Data=Data,simplify=FALSE)
  TestA= sapply(flds,ArrayTest,Data=Data,simplify=FALSE)
  SUB=sapply(flds,subDist,D=D,simplify=FALSE)
  hc=sapply(SUB,mycluster2,k=k)
  C=mapply(kcentroid,TrainA,hc,SIMPLIFY=FALSE)
  CentrCor= mapply(cor,C[1],C[2:nCV],SIMPLIFY = FALSE)
  Cord= sapply(CentrCor,alloc2,simplify=FALSE)
  Cord=append(list(as.matrix(1:k,ncol=1)),Cord , 1)
  C=mapply(reord,C,Cord,SIMPLIFY=FALSE)
  CLASS=mapply(classify2,TestA,C,SIMPLIFY=FALSE)
  #t1=Sys.time();for(i in 1:10){CLASS=mapply(classify2,TestA,C,SIMPLIFY=FALSE)};t2=Sys.time();t2-t1
  clustclass=unlist(CLASS)
  clustclass=clustclass[order(as.vector(unlist(flds)))]
  fit1=mean(silhouette(clustclass,D)[,3])
  fit2= fitness(OS,clustclass)
  
  return(c(fit1,fit2))
}

#Minimum number of genes to use in a solution (A constraint for the algorithm)
minGenes= function(x){ sum(x)>=10 & sum(x)< chrom_length}

#http://ictactjournals.in/paper/IJSC_V6_I1_paper_4_pp_1083_1092.pdf
#Multiple point crossover
#a and b is solution 1 and 2 respectively (binary vectors)
#n is the number of cut points
kcrossover=function(a,b,n){
  if(length(a)!=length(b)) {stop("vectors of unequal length")}
  l=length(a)
  if(n>=(length(a)-1)) {stop("number of cut points bigger than possible sites")}
  points=sample(2:(l-1),n,replace=FALSE)
  to=c(points[order(points)][-n],l)
  from= c(1,to[-length(to)]+1)
  cutpoints=list()
  for(i in 1:n){
    cutpoints[[i]]=seq(from[i],to[i])
  }
  achild=as.numeric()
  bchild=as.numeric()
  for(i in 1:n){
    if(i%%2==0){
      achild=c(achild,a[cutpoints[[i]]])
      bchild=c(bchild,b[cutpoints[[i]]])
    }else{
      achild=c(achild,b[cutpoints[[i]]])
      bchild=c(bchild,a[cutpoints[[i]]])
    }}
  return(list(achild,bchild))
}

#Uniformcrossover
ucrossover=function(a,b){
  if(length(a)!=length(b)) {stop("vectors of unequal length")}
  l=length(a)
  points=as.logical(rbinom(l,1,prob=runif(1)))
  achild=numeric(l)
  achild[points]=a[points]
  achild[!points]=b[!points]
  bchild=numeric(l)
  bchild[points]=b[points]
  bchild[!points]=a[!points]
  return(list(achild,bchild))
}

#Asymmetric mutation operator:
#Analysis of an Asymmetric Mutation Operator; Jansen et al.

asMut=function(x){
  chrom_length=length(x)
  res=x
  Active=sum(x)
  Deactive=chrom_length-Active
  mutrate1=1/Active
  mutrate2=1/Deactive
  mutpoint1=sample(c(1,0),Deactive,prob=c(mutrate2,1-mutrate2),replace=TRUE)
  res[x==0]=abs(x[x==0]-mutpoint1)
  
  mutpoint2=sample(c(1,0),Active,prob=c(mutrate1,1-mutrate1),replace=TRUE)
  res[x==1]=abs(x[x==1]-mutpoint2)
  
  return(res)
}


#Offspring creation
offspring= function(X1=X, Nclust=Nclust){
  
  New= matrix(NA,ncol=chrom_length,nrow=population) #Create empty matrix to add new individuals
  NewK= matrix(NA,nrow=1,ncol=population)           #same for cluster chromosome
  
  matingPool <- tournamentSelection(X1,population,TournamentSize) #Use tournament selection, to select parents that will give offspring
  
  count=0 #Count how many offsprings are still needed to reach the original population size
  while(anyNA(New)){
    count=count+1
    ##a.Select a pair of parent chromosomes from the matingPool
    Pair=sample(1: nrow(matingPool),2,replace=F)
    
    
    ##b.With probability pc (the "crossover probability" or "crossover rate"), cross over the pair at a n randomly chosen points (with probability p, chosen randomly from uniform distribution) to form two offspring. If no crossover takes place, exact copies of their respective parents are pass to the next generation.
    
    Cp=1    # with elitism there is no need to add a crossover probability 
    if(sample(c(1,0),1,p=c(Cp,1-Cp))== 1){
      #multiple point crossingover    
      offspring=ucrossover(matingPool[Pair[1],1:chrom_length],matingPool[Pair[2],1:chrom_length])
      off1= offspring[[1]]
      off2= offspring[[2]]
    }else{
      off1= matingPool[Pair[1],1:chrom_length]
      off2= matingPool[Pair[2],1:chrom_length]
    }
    
    ##c.Mutate the two offspring at each locus with probability Mp (the mutation probability or mutation rate),
    # and place the resulting chromosomes in the new population.
    #since the results are sparse strings, cosine similarity is more adequate
    #Mutation by asymmetric mutation: Analysis of an Asymmetric Mutation Operator; Jansen et al.
    
    Mp=cosine(matingPool[Pair[1],1:chrom_length],matingPool[Pair[2],1:chrom_length])
    
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      off1=asMut(off1)
    }
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      off2=asMut(off2)
    }
    
    #Mutation for k (number of partitions)
    
    p=0.3
    probs= rep((1-p)/9,9)
    k1=matingPool[Pair[1],"k"]
    k2=matingPool[Pair[2],"k"]
    probs[k1-1]=probs[k1-1]+p/2;probs[k2-1]=probs[k2-1]+p/2
    offk1= sample(2:10,1,prob=probs)
    offk2= sample(2:10,1,prob=probs)

    #Add offsprings to new generation
    New[count,]= off1
    New[count+1,]= off2
    NewK[,count]=offk1
    NewK[,count+1]=offk2
  } 
  return(list(New=New,NewK=NewK))
}

pen=function(x){ 1 / (1+(x/500)^2)}

#Stable allocation of clusters (gale shapley alg)
alloc=function(C1,C2){
  x1=as.matrix(gpuCor(C1,C2))
  x2=t(x1)
  ord=galeShapley.marriageMarket(x1,x2)
  return(ord$engagements)
}
