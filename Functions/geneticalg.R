
#parallele computing
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

#Option to run it in GPU mode or CPU mode
if(GPU==TRUE){
  gpupackage="gpuR"
  computingtype= "using GPU computing"
}else{
  gpupackage="amap"
  mydist<<-function(c) {Dist(c,method="correlation")}
  computingtype="using CPU computing"
}

PARETO=list() #Empty list to save the solutions

#1. Create random population of solutions

#Creating random clusters from 2-10
Nclust= sample(2:10,population,replace=TRUE)

#Matrix with random TRUE false with uniform distribution, representing solutions to test
X=matrix(NA,nrow=population,ncol=chrom_length)
for(i in 1:population){
  prob=runif(1,0,1)
  X[i,]= sample(c(1,0),chrom_length,replace=T,prob=c(prob,1-prob))
}

#####Main loop
for(g in 1:generations){
  
  start_time <- Sys.time() #Measures generation time

#2.Calculate the fitness f(x) of each chromosome x in the population.
  Fit1= apply(X,1,minGenes) #Apply constraints (min 10 genes per solution)
  X= X[Fit1,] 
  X=apply(X,2,as.logical)
  n=nrow(X)
  Nclust=Nclust[Fit1]
  
  k=Nclust
  
  flds = createFolds(1:ncol(prob_matrix),k=nCV)
  
  #Calculate Fitnes 1 (silhouette) and 2 (Survival differences)
   Fit2=foreach(i=1:nrow(X),.packages=c('cluster',gpupackage,"cba","survival","matchingR"),.combine=rbind) %dopar% {
    crossvalidation(prob_matrix,flds,X[i,],k[i],OS=OS)
  }   
  
  # penalization of SC by number of genes
  Fit2[,1]=(Fit2[,1]*pen(rowSums(X))) 
  
  if(g==1){
    
    PARETO[[g]]= Fit2 #Saves the fitnes of the solutions of the current generation 
    ranking=fastNonDominatedSorting(Fit2*-1) #NonDominatedSorting from package nsga2R. -1 multiplication to alter the sign of the fitness (function sorts from min to max)
    rnkIndex <- integer(n)
    i <- 1
    while (i <= length(ranking)) { #saves the rank of each solution
      rnkIndex[ranking[[i]]] <- i
      i <- i + 1
    }
    
    X1=cbind(X,k,Fit2,rnkIndex) #data.frame with solution vector, number of clusters and ranking
    objRange=apply(Fit2,2,max)-apply(Fit2,2,min) # Range of fitness of the solutions 
    CrowD= crowdingDist4frnt(X1,ranking,objRange) #Crowding distance of each front (nsga2R package)
    CrowD=apply(CrowD,1,sum)
    
    X1= cbind(X1,CrowD) #data.frame with solution vector, number of clusters, ranking and crowding distance
    
    Archive[[g]]=X1[X1[,"rnkIndex"]==1,(chrom_length+2):(chrom_length+3)]
    #output for the generation
    print(paste0("Generation ",g," Non-dominated solutions:"))  
    print(X1[X1[,"rnkIndex"]==1,(chrom_length+1):(chrom_length+5)])
    
    #save parent generation
    Xold=X
    Nclustold=Nclust
    
    #3. create offspring
    NEW= offspring(X1)
    X=NEW[["New"]]
    Nclust=NEW[["NewK"]]
    
  }else{
    
    oldnew=rbind(PARETO[[g-1]],Fit2)
    oldnewfeature= rbind(Xold,X)
    oldnewNclust= c(Nclustold,Nclust)
    oldnewK=oldnewNclust
    
    ranking2=fastNonDominatedSorting(oldnew*-1) #NonDominatedSorting from package nsga2R. -1 multiplication to alter the sign of the fitness (function sorts from min to max)
    rnkIndex2 <- integer(nrow(oldnew))
    i <- 1
    while (i <= length(ranking2)) { #saves the rank of each solution
      rnkIndex2[ranking2[[i]]] <- i
      i <- i + 1
    }
    
    
    X1=cbind(oldnewfeature,k=oldnewK,oldnew,rnkIndex=rnkIndex2)
    objRange=apply(oldnew,2,max)-apply(oldnew,2,min) # Range of fitness of the solutions 
    CrowD= crowdingDist4frnt(X1,ranking2,objRange) #Crowding distance of each front (nsga2R package)
    CrowD=apply(CrowD,1,sum)
    X1= cbind(X1,CrowD) #data.frame with solution vector, number of clusters, ranking and crowding distance
    X1= X1[X1[,"CrowD"]>0,]
    O=order(X1[,"rnkIndex"],X1[,"CrowD"]*-1)[1:population]
    X1=X1[O,]
    
    PARETO[[g]]= X1[,(chrom_length+2):(chrom_length+3)] #Saves the fitnes of the solutions of the current generation 
    
    #output for the generation
    print(paste0("Generation ",g," Non-dominated solutions:"))  
    print(X1[X1[,"rnkIndex"]==1,(chrom_length+1):(chrom_length+5)])
    
    Xold=X1[,1:chrom_length]
    Nclustold=X1[,"k"]
    
    NEW= offspring(X1)
    X=NEW[["New"]]
    Nclust=NEW[["NewK"]]
  }
  
  #5.Go to step 2
  end_time <- Sys.time()
  t=end_time-start_time
  
  print(t)
  gc()
  
  if(g%%50==0){
    colnames(X1)[1:(ncol(X1)-5)]=rownames(prob_matrix)
    output=list(Solutions=X1,ParetoFront=PARETO)
    filename= paste0(resultdir,"/results",g,".rda")
    save(file=filename,output)
    }
  
}

stopCluster(cluster)

