##################
####Functions#####
##################

#Distance type used in the algorithm

distancetype="pearson"          #Options are: "pearson","uncentered" (cosine),"spearman","euclidean"


#Clinical Tidy METABRIC and OSLO data

NotInOslo_NAimput= function(clinic){
  idx = grep("NOT_IN_OSLOVAL", colnames(clinic))
  if (length(idx) > 0) 
    clinic = clinic[, -idx]
  idx = grep("IHC", colnames(clinic))
  if (length(idx) > 0) 
    clinic = clinic[, -idx]
  for (i in 1:ncol(clinic)) {
    idx = which(is.na(clinic[, i]))
    if (length(idx) > 0) {
      if (class(clinic[, i]) == "numeric") 
        clinic[idx, i] = mean(clinic[-idx, i])
      if (class(clinic[, i]) == "factor") {
        cc = as.vector(clinic[, i])
        cc[idx] = "NA"
        clinic[, i] = factor(cc)
      }
    }
  }
  return(clinic)
}

reshapeClinic= function(clinic)  {
  h.IDC = as.numeric(clinic$histological_type == "IDC")
  h.ILC = as.numeric(clinic$histological_type == "ILC")
  h.IDCpILC = as.numeric(clinic$histological_type == "IDC+ILC")
  h.IDCnMED = as.numeric(clinic$histological_type == "IDC-MED")#Not 1 in Oslo
  h.IDCnMUC = as.numeric(clinic$histological_type == "IDC-MUC")
  h.IDCnTUB = as.numeric(clinic$histological_type == "IDC-TUB")
  er.P = as.numeric(clinic$ER.Expr == "+")
  er.N = as.numeric(clinic$ER.Expr == "-")
  tr.CT = as.numeric((clinic$Treatment == "CT") | (clinic$Treatment == "CT/HT") | (clinic$Treatment == "CT/HT/RT") | (clinic$Treatment == "CT/RT"))
  tr.HT = as.numeric((clinic$Treatment == "HT") | (clinic$Treatment == "CT/HT") | (clinic$Treatment == "HT/RT") | (clinic$Treatment == "CT/HT/RT"))
  tr.RT = as.numeric((clinic$Treatment == "RT") | (clinic$Treatment == "CT/HT/RT") | (clinic$Treatment == "CT/RT") | (clinic$Treatment == "HT/RT"))
  gd.1 = as.numeric(clinic$grade == 1)
  gd.2 = as.numeric(clinic$grade == 2)
  gd.3 = as.numeric(clinic$grade == 3)
  cmat <- data.frame(clinic[, c(1:3)], gd.1, gd.2, gd.3, h.IDC, h.ILC, h.IDCpILC, h.IDCnMED, h.IDCnMUC, h.IDCnTUB, er.N,er.P, tr.CT, tr.HT, tr.RT)
  for (i in 4:ncol(cmat)) {
    cmat[, i] = factor(cmat[, i])
  }
  rownames(cmat)=make.names(rownames(cmat))
  return(cmat)
}


#Expression set tidy for METABRIC and OSLO
expressiontidy2= function(Expression,ANNOT, CLINIC, SURV,BATCH){
  #require(illuminaio)
  require(AnnotationDbi)
  #require(limma)
  
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



distancetype=match.arg(distancetype,c("pearson","uncentered","euclidean"))
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

#There are papers saing this is the best for gene
#https://core.ac.uk/download/pdf/61320037.pdf
#Also usefull for cosine distance (is the same)
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
#1/9/18 cambié parametros do.swap=FALSE to do.swap=TRUE
mycluster2<- function(c,k){
  return(list(cluster=pam(c,k,cluster.only=TRUE,diss=TRUE,do.swap=TRUE,keep.diss=FALSE,keep.data=FALSE,pamonce=2)))
}

#Humming distance

humdist= function(x,y){sum(abs(x-y))}   #Function to calculate Humming distance between two binary vectors

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


#classify2=function(data,centroid){
#  R=cor(data,centroid,method=distancetype)
#  scores<-apply(R,1,which.max)
#  return(scores)
#}

#crossvalidation function based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3105299/
#Data is the expression matrix
#flds is a list with the indexes to partition the data in n different folds
#indv is the solution to test, namely, a binary vector to subset the genes to test
#vit_stat is the codification for an event in the clinical data (e.g. "Dead"), is used in the fitness function

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

#Modified Binary lector to read cluster codification (3 bit binary)
BinToDec <- function(x) 
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))+2

DecToBin<- function(x){
  x=x-2
  rev(as.integer(intToBits(x))[1:3])
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
    
    
    ##b.With probability pc (the "crossover probability" or "crossover rate"), cross over the pair at a
    #n randomly chosen points (with probability p, chosen randomly from uniform distribution) to form two offspring. If no
    #crossover takes place, exact copies of their respective parents are pass to the next generation.
    
    Cp=1    # Con el esquema del elitismo no es necesario agregarle probabilidad de no hacer crossover
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
    #HD= humdist(matingPool[Pair[1],1:chrom_length],matingPool[Pair[2],1:chrom_length])
    #Mp= (chrom_length-HD)/chrom_length
    #since the results are sparse strings, cosine similarity is more adequate
    #Mutation by asymmetric mutation: Analysis of an Asymmetric Mutation Operator; Jansen et al.
    
    Mp=cosine(matingPool[Pair[1],1:chrom_length],matingPool[Pair[2],1:chrom_length])
    
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      off1=asMut(off1)
    }
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      off2=asMut(off2)
    }
    
    #Mutation for cluster N
    
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

subfeatures= function(prob_matrix,Nfeatures,Sigvector,genesig){
  rgenes= Nfeatures-length(Sigvector)
  nonsiggenes= genesig[-match(Sigvector, genesig)]
  set.seed(412)
  randomgenes=sample(nonsiggenes,rgenes,replace=FALSE)
  finalsig= c(as.character(Sigvector),as.character(randomgenes))
  return(prob_matrix[finalsig,])
}



pen=function(x){ 1 / (1+(x/500)^2)}

#Stable allocation of clusters (gale shapley alg)
alloc=function(C1,C2){
  x1=as.matrix(gpuCor(C1,C2))
  x2=t(x1)
  ord=galeShapley.marriageMarket(x1,x2)
  return(ord$engagements)
}

#Heatmap.3 function

heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
        
        invalid <- function (x) {
                if (missing(x) || is.null(x) || length(x) == 0)
                        return(TRUE)
                if (is.list(x))
                        return(all(sapply(x, invalid)))
                else if (is.vector(x))
                        return(all(is.na(x)))
                else return(FALSE)
        }
        
        x <- as.matrix(x)
        scale01 <- function(x, low = min(x), high = max(x)) {
                x <- (x - low)/(high - low)
                x
        }
        retval <- list()
        scale <- if (symm && missing(scale))
                "none"
        else match.arg(scale)
        dendrogram <- match.arg(dendrogram)
        trace <- match.arg(trace)
        density.info <- match.arg(density.info)
        if (length(col) == 1 && is.character(col))
                col <- get(col, mode = "function")
        if (!missing(breaks) && (scale != "none"))
                warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                        "specified can produce unpredictable results.", "Please consider using only one or the other.")
        if (is.null(Rowv) || is.na(Rowv))
                Rowv <- FALSE
        if (is.null(Colv) || is.na(Colv))
                Colv <- FALSE
        else if (Colv == "Rowv" && !isTRUE(Rowv))
                Colv <- FALSE
        if (length(di <- dim(x)) != 2 || !is.numeric(x))
                stop("`x' must be a numeric matrix")
        nr <- di[1]
        nc <- di[2]
        if (nr <= 1 || nc <= 1)
                stop("`x' must have at least 2 rows and 2 columns")
        if (!is.numeric(margins) || length(margins) != 2)
                stop("`margins' must be a numeric vector of length 2")
        if (missing(cellnote))
                cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
        if (!inherits(Rowv, "dendrogram")) {
                if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                             c("both", "row"))) {
                        if (is.logical(Colv) && (Colv))
                                dendrogram <- "column"
                        else dedrogram <- "none"
                        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                                dendrogram, "'. Omitting row dendogram.")
                }
        }
        if (!inherits(Colv, "dendrogram")) {
                if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                             c("both", "column"))) {
                        if (is.logical(Rowv) && (Rowv))
                                dendrogram <- "row"
                        else dendrogram <- "none"
                        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                                dendrogram, "'. Omitting column dendogram.")
                }
        }
        if (inherits(Rowv, "dendrogram")) {
                ddr <- Rowv
                rowInd <- order.dendrogram(ddr)
        }
        else if (is.integer(Rowv)) {
                hcr <- hclustfun(distfun(x))
                ddr <- as.dendrogram(hcr)
                ddr <- reorder(ddr, Rowv)
                rowInd <- order.dendrogram(ddr)
                if (nr != length(rowInd))
                        stop("row dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Rowv)) {
                Rowv <- rowMeans(x, na.rm = na.rm)
                hcr <- hclustfun(distfun(x))
                ddr <- as.dendrogram(hcr)
                ddr <- reorder(ddr, Rowv)
                rowInd <- order.dendrogram(ddr)
                if (nr != length(rowInd))
                        stop("row dendrogram ordering gave index of wrong length")
        }
        else {
                rowInd <- nr:1
        }
        if (inherits(Colv, "dendrogram")) {
                ddc <- Colv
                colInd <- order.dendrogram(ddc)
        }
        else if (identical(Colv, "Rowv")) {
                if (nr != nc)
                        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
                if (exists("ddr")) {
                        ddc <- ddr
                        colInd <- order.dendrogram(ddc)
                }
                else colInd <- rowInd
        }
        else if (is.integer(Colv)) {
                hcc <- hclustfun(distfun(if (symm)
                        x
                        else t(x)))
                ddc <- as.dendrogram(hcc)
                ddc <- reorder(ddc, Colv)
                colInd <- order.dendrogram(ddc)
                if (nc != length(colInd))
                        stop("column dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Colv)) {
                Colv <- colMeans(x, na.rm = na.rm)
                hcc <- hclustfun(distfun(if (symm)
                        x
                        else t(x)))
                ddc <- as.dendrogram(hcc)
                ddc <- reorder(ddc, Colv)
                colInd <- order.dendrogram(ddc)
                if (nc != length(colInd))
                        stop("column dendrogram ordering gave index of wrong length")
        }
        else {
                colInd <- 1:nc
        }
        retval$rowInd <- rowInd
        retval$colInd <- colInd
        retval$call <- match.call()
        x <- x[rowInd, colInd]
        x.unscaled <- x
        cellnote <- cellnote[rowInd, colInd]
        if (is.null(labRow))
                labRow <- if (is.null(rownames(x)))
                        (1:nr)[rowInd]
        else rownames(x)
        else labRow <- labRow[rowInd]
        if (is.null(labCol))
                labCol <- if (is.null(colnames(x)))
                        (1:nc)[colInd]
        else colnames(x)
        else labCol <- labCol[colInd]
        if (scale == "row") {
                retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
                x <- sweep(x, 1, rm)
                retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
                x <- sweep(x, 1, sx, "/")
        }
        else if (scale == "column") {
                retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
                x <- sweep(x, 2, rm)
                retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
                x <- sweep(x, 2, sx, "/")
        }
        if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
                if (missing(col) || is.function(col))
                        breaks <- 16
                else breaks <- length(col) + 1
        }
        if (length(breaks) == 1) {
                if (!symbreaks)
                        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                                      length = breaks)
                else {
                        extreme <- max(abs(x), na.rm = TRUE)
                        breaks <- seq(-extreme, extreme, length = breaks)
                }
        }
        nbr <- length(breaks)
        ncol <- length(breaks) - 1
        if (class(col) == "function")
                col <- col(ncol)
        min.breaks <- min(breaks)
        max.breaks <- max(breaks)
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks
        if (missing(lhei) || is.null(lhei))
                lhei <- c(keysize, 4)
        if (missing(lwid) || is.null(lwid))
                lwid <- c(keysize, 4)
        if (missing(lmat) || is.null(lmat)) {
                lmat <- rbind(4:3, 2:1)
                
                if (!missing(ColSideColors)) {
                        #if (!is.matrix(ColSideColors))
                        #stop("'ColSideColors' must be a matrix")
                        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                                stop("'ColSideColors' must be a matrix of nrow(x) rows")
                        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
                        #lhei <- c(lhei[1], 0.2, lhei[2])
                        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
                }
                
                if (!missing(RowSideColors)) {
                        #if (!is.matrix(RowSideColors))
                        #stop("'RowSideColors' must be a matrix")
                        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                                stop("'RowSideColors' must be a matrix of ncol(x) columns")
                        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
                        #lwid <- c(lwid[1], 0.2, lwid[2])
                        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
                }
                lmat[is.na(lmat)] <- 0
        }
        
        if (length(lhei) != nrow(lmat))
                stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
        if (length(lwid) != ncol(lmat))
                stop("lwid must have length = ncol(lmat) =", ncol(lmat))
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        
        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
        
        if (!missing(RowSideColors)) {
                if (!is.matrix(RowSideColors)){
                        par(mar = c(margins[1], 0, 0, 0.5))
                        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
                } else {
                        par(mar = c(margins[1], 0, 0, 0.5))
                        rsc = t(RowSideColors[,rowInd, drop=F])
                        rsc.colors = matrix()
                        rsc.names = names(table(rsc))
                        rsc.i = 1
                        for (rsc.name in rsc.names) {
                                rsc.colors[rsc.i] = rsc.name
                                rsc[rsc == rsc.name] = rsc.i
                                rsc.i = rsc.i + 1
                        }
                        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
                        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
                        if (length(rownames(RowSideColors)) > 0) {
                                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
                        }
                }
        }
        
        if (!missing(ColSideColors)) {
                
                if (!is.matrix(ColSideColors)){
                        par(mar = c(0.5, 0, 0, margins[2]))
                        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
                } else {
                        par(mar = c(0.5, 0, 0, margins[2]))
                        csc = ColSideColors[colInd, , drop=F]
                        csc.colors = matrix()
                        csc.names = names(table(csc))
                        csc.i = 1
                        for (csc.name in csc.names) {
                                csc.colors[csc.i] = csc.name
                                csc[csc == csc.name] = csc.i
                                csc.i = csc.i + 1
                        }
                        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
                        image(csc, col = as.vector(csc.colors), axes = FALSE)
                        if (length(colnames(ColSideColors)) > 0) {
                                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
                        }
                }
        }
        
        par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)
        cellnote <- t(cellnote)
        if (revC) {
                iy <- nr:1
                if (exists("ddr"))
                        ddr <- rev(ddr)
                x <- x[, iy]
                cellnote <- cellnote[, iy]
        }
        else iy <- 1:nr
        image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
        retval$carpet <- x
        if (exists("ddr"))
                retval$rowDendrogram <- ddr
        if (exists("ddc"))
                retval$colDendrogram <- ddc
        retval$breaks <- breaks
        retval$col <- col
        if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
                mmat <- ifelse(is.na(x), 1, NA)
                image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
                      col = na.color, add = TRUE)
        }
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
             cex.axis = cexCol)
        if (!is.null(xlab))
                mtext(xlab, side = 1, line = margins[1] - 1.25)
        axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
             cex.axis = cexRow)
        if (!is.null(ylab))
                mtext(ylab, side = 4, line = margins[2] - 1.25)
        if (!missing(add.expr))
                eval(substitute(add.expr))
        if (!missing(colsep))
                for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
        if (!missing(rowsep))
                for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
        min.scale <- min(breaks)
        max.scale <- max(breaks)
        x.scaled <- scale01(t(x), min.scale, max.scale)
        if (trace %in% c("both", "column")) {
                retval$vline <- vline
                vline.vals <- scale01(vline, min.scale, max.scale)
                for (i in colInd) {
                        if (!is.null(vline)) {
                                abline(v = i - 0.5 + vline.vals, col = linecol,
                                       lty = 2)
                        }
                        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
                        xv <- c(xv[1], xv)
                        yv <- 1:length(xv) - 0.5
                        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
                }
        }
        if (trace %in% c("both", "row")) {
                retval$hline <- hline
                hline.vals <- scale01(hline, min.scale, max.scale)
                for (i in rowInd) {
                        if (!is.null(hline)) {
                                abline(h = i + hline, col = linecol, lty = 2)
                        }
                        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
                        yv <- rev(c(yv[1], yv))
                        xv <- length(yv):1 - 0.5
                        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
                }
        }
        if (!missing(cellnote))
                text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
                     col = notecol, cex = notecex)
        par(mar = c(margins[1], 0, 0, 0))
        if (dendrogram %in% c("both", "row")) {
                plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
        }
        else plot.new()
        par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
        if (dendrogram %in% c("both", "column")) {
                plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
        }
        else plot.new()
        if (!is.null(main))
                title(main, cex.main = 1.5 * op[["cex.main"]])
        if (key) {
                par(mar = c(5, 4, 2, 1), cex = 0.75)
                tmpbreaks <- breaks
                if (symkey) {
                        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
                        min.raw <- -max.raw
                        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
                        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
                }
                else {
                        min.raw <- min(x, na.rm = TRUE)
                        max.raw <- max(x, na.rm = TRUE)
                }
                
                z <- seq(min.raw, max.raw, length = length(col))
                image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                      xaxt = "n", yaxt = "n")
                par(usr = c(0, 1, 0, 1))
                lv <- pretty(breaks)
                xv <- scale01(as.numeric(lv), min.raw, max.raw)
                axis(1, at = xv, labels = lv)
                if (scale == "row")
                        mtext(side = 1, "Row Z-Score", line = 2)
                else if (scale == "column")
                        mtext(side = 1, "Column Z-Score", line = 2)
                else mtext(side = 1, KeyValueName, line = 2)
                if (density.info == "density") {
                        dens <- density(x, adjust = densadj, na.rm = TRUE)
                        omit <- dens$x < min(breaks) | dens$x > max(breaks)
                        dens$x <- dens$x[-omit]
                        dens$y <- dens$y[-omit]
                        dens$x <- scale01(dens$x, min.raw, max.raw)
                        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                              lwd = 1)
                        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
                        title("Color Key\nand Density Plot")
                        par(cex = 0.5)
                        mtext(side = 2, "Density", line = 2)
                }
                else if (density.info == "histogram") {
                        h <- hist(x, plot = FALSE, breaks = breaks)
                        hx <- scale01(breaks, min.raw, max.raw)
                        hy <- c(h$counts, h$counts[length(h$counts)])
                        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                              col = denscol)
                        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
                        title("Color Key\nand Histogram")
                        par(cex = 0.5)
                        mtext(side = 2, "Count", line = 2)
                }
                else title("Color Key")
        }
        else plot.new()
        retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                        high = retval$breaks[-1], color = retval$col)
        invisible(retval)
}
