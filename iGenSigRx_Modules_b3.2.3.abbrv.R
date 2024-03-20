require(stringr)
require(vegan)
require(qvalue)
require(tidyr)
require(moments)
require(limma)
require(dplyr)
require(gdata)
require(survival)
require(pROC)
require(preprocessCore)
require(DescTools)
require(ggplot2)
require(robustbase)
require(pheatmap)
require(readr)
require(fastcluster)
require(Rfast)
require(dynamicTreeCut)
require(data.table)

##########################################################
#########       General Modules                ###########
##########################################################

#scale the max of a vector to 1 and min to 0
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

gm_mean = function(x, na.rm=TRUE){
  if (sum(x)==0){
    return(0)
  }else{
    #return(exp(mean(log(x))))
    return(exp(sum(log(x[x>0]), na.rm=na.rm)/length(x))) 
  }
}

#filter by one list
filter <- function (x,filter,inlist) {
  if (inlist){
    x <- x[x %in% filter]  
  }
  else{
    x <- x[!x %in% filter]  
  }
  return (x)
}


read_gmt <- function(fileName, min = 5) {
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  read.list <- lapply (Lines, function(x) {
    genotype=strsplit(x, "\t")[[1]]
    gname=paste(genotype[1],genotype[2],sep=":")
    return(list(gname,genotype[-c(1:2)]))
  })
  genotype.list=lapply(read.list, `[[`, 2) 
  names(genotype.list)= lapply(read.list, `[[`, 1)
  genotype.list=genotype.list[lengths(genotype.list)>=min]
  return(genotype.list)
}

read_preCal <- function(fileName) {
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  preCal.list=list()
  for (Line in Lines){
    tmp.line<-unlist(strsplit(Line,"\t"))
    tmp.split=strsplit(tmp.line[3:length(tmp.line)],"@")
    feature.redun=as.numeric(sapply(tmp.split, `[[`, 2))
    names(feature.redun)= sapply(tmp.split, `[[`, 1)
    preCal.list[[tmp.line[1]]]=feature.redun
  }
  return(preCal.list)
}

#find concepts that contains the cosmicID in a conceptList
findConcepts<-function(cosmicID,conceptList){
  concepts<-c()
  for(i in 1:length(conceptList)){
    if(cosmicID %in% conceptList[[i]]){
      concepts<-append(concepts,names(conceptList)[i])
    }
  }
  return(concepts)
}

#turn genotype list into a matrix fast
genolist2matrix<-function(genolist){
  un.genolist <- unlist(genolist)
  cell<-sort(unique(un.genolist))
  tmp.bin <- matrix(0, nrow = length(genolist), ncol = length(cell))
  dimnames(tmp.bin) <- list(names(genolist), cell)
  ij <- cbind(rep(1:length(genolist), lengths(genolist)), match(un.genolist,cell))
  tmp.bin[ij] <- 1
  return(tmp.bin)
}

#turn matrix into genolist file
matrix2genofile<-function(matrix,outfile){
  for (i in 1:nrow(matrix)){
    tmp.out=c(unlist(strsplit(rownames(matrix)[i],":")),names(which(matrix[i,]==1)))
    write(tmp.out,file=outfile,append=TRUE,sep="\t",ncolumns = length(tmp.out))
  }
}
preCalList2gmt<-function(preCal.list){
  preCal.matrix=matrix(ncol = 1, nrow = length(preCal.list))
  for (i in 1:length(preCal.list)){
    ECN=length(preCal.list[[i]])/mean(as.numeric(preCal.list[[i]]),trim=0.3)
    redundancy=paste(paste(names(preCal.list[[i]]),preCal.list[[i]],sep="@"),collapse="\t")
    preCal.matrix[i,]=paste(names(preCal.list)[i],ECN,redundancy,sep="\t")
  }
  return(preCal.matrix)
}
#simulate sequencing errors in gmt files
simulate.seqError<-function(genotypefile,outpath,errPerc=0.05,permuN=5){
  genotype.list<-read_gmt(genotypefile,min=5)
  genotype.mat<-genolist2matrix(genotype.list)
  for (permu in 1:permuN){
    tmp.mat=apply(genotype.mat,2,function(x){
      index=sample(1:nrow(genotype.mat), round(errPerc * nrow(genotype.mat)))
      x[index]=1-x[index]
      return(x)
    })
    rownames(tmp.mat)=rownames(genotype.mat)
    permu.diff=sum(tmp.mat!=genotype.mat)/length(genotype.mat)
    print (paste("Permutation:",permu,"differences:",permu.diff))
    outfile=paste(outpath,"/",sub(pattern = ".gmt", replacement = "", basename(genotypefile)),".permu",permu,".gmt",sep="")
    writeLines(paste("#Simulating",errPerc*100,"percent sequencing error rates",sep=" "),outfile)
    matrix2genofile(matrix=tmp.mat,outfile = outfile)
  }
}
##########################################################
#########   Preparation Modules                ###########
##########################################################

#create sample training (0) and testing sets (1) based on n folds
FoldCrossVal.dataset<-function(sample,fold=5){
    sample <- sample[sample(length(sample))]
    folds <- cut(seq(1,length(sample)),breaks=fold,labels=FALSE)
    sample.assign=data.frame()
    for(i in 1:fold){
      #Segement your data by fold using the which() function 
      testSet <- sample[which(folds==i,arr.ind=TRUE)]
      sample.assign=rbind(sample.assign, ifelse(sample %in% testSet,1,0))
    }
    colnames(sample.assign)=sample
    return(sample.assign)
}
#create sample training (0) and testing sets (1) based on random sampling
SampleCrossVal.dataset<-function(phenoData,sample.col,label.col,test.perc,repeatn=5,confound=NULL){ #label.col: column for categorical label
  sample.assign=data.frame()
  label.levels=unique(phenoData[,label.col])
  for(i in 1:repeatn){  
    subject.test=c()
    for (level in label.levels){
      if (typeof(confound)=="NULL"){
        level.pos=which(phenoData[,label.col]==level)
        subject.test=c(subject.test,phenoData[,sample.col][sample(level.pos,round(test.perc*length(level.pos)))]) 
      }else{
        cf.levels=unique(phenoData[,confound])
        for (cf.level in cf.levels){
          cf.pos=which(phenoData[,confound]==cf.level & phenoData[,label.col]==level)
          subject.test=c(subject.test,phenoData[,sample.col][sample(cf.pos,round(test.perc*length(cf.pos)))]) 
          }
      }
    }
   sample.assign=rbind(sample.assign, ifelse(phenoData[,sample.col] %in% subject.test,1,0))
  }
  colnames(sample.assign)=phenoData[,sample.col]
  return(sample.assign)
}
#create sample training (0) and testing sets (1) based on random sampling and equal responder/nonresponder
SampleCrossVal.response<-function(phenoData,sample.col,label.col,responder.train,repeatn=5){ #label.col: column for categorical label
  if (sum(phenoData[,label.col] %in% c("sen","res"))!=nrow(phenoData)){
    die("label column must contain 'sen' or 'res' values only")
  }
  sample.assign=data.frame()
  label.levels=unique(phenoData[,label.col])
  for(i in 1:repeatn){  
    subject.test=c()
    train.sen=phenoData[,sample.col][sample(which(phenoData[,label.col]=="sen"),responder.train)]
    train.res=phenoData[,sample.col][sample(which(phenoData[,label.col]=="res"),responder.train)]
    sample.assign=rbind(sample.assign, ifelse(phenoData[,sample.col] %in% c(train.sen,train.res),0,1))
  }
  colnames(sample.assign)=phenoData[,sample.col]
  return(sample.assign)
}

#Freeze version (Trim Corrected): Generate 6 Overlapping Expression Bins Generated Based On Cutoffs created using trimmed mean and sd of fold change using unlog TPM expression data#######################################
#sdlevels means the levels of standard deviations must be symmetric and must include zero (the default is -3 to 3)
binarize.expfeature<-function(expData,outfilepath,sdpercfilter=0.2,sdlevels=-6:6,genomewide=TRUE){
  suppressWarnings(file.remove(outfilepath))
  if (length(grep("_1$", colnames(expData)))>0){
    expData<-expData[, -grep("_1$", colnames(expData))]#remove duplcate columns (Duplicated column names deduplicated: '1503362' => '1503362_1' [845])  
  }
  rownames(expData)<-gsub("\\.\\d+$", "", rownames(expData)) #remove .\d+ in ensemble_gene
  expMatrix<-as.matrix(expData)
  if (genomewide==TRUE){
    expMatrix=normalize.quantiles(expMatrix)
    dimnames(expMatrix) <- list(rownames(expData), colnames(expData))
    exp.sd <- apply(expMatrix, 1, sd)
    expMatrix<-expMatrix[which(exp.sd > quantile(exp.sd, sdpercfilter)), ] 
  }
  samples=colnames(expMatrix)
  writeLines(paste(c(paste("#Overlapping Expression Bins Generated Based On Cutoffs created using Trimmed (0.1) Mean + ",paste(sdlevels[sdlevels!=0],collapse="/"),"*SD based on log2 transformed fold change for quantile normalized TPM",sep=""),samples),collapse='\t'),outfilepath)
  for (i in 1:nrow(expMatrix)){
    expi=as.numeric(expMatrix[i,])
    names(expi)=names(expMatrix[i,])
    #expi=expi+min(expi[expi!=0])/2
    expi=expi+1
    expi=log2(expi/mean(expi,trim=0.1))
    expi.trim=Trim(expi,trim=0.1)
    if (sd(expi.trim)==0){
      next
    }
    cutoffs=mean(expi.trim)+sd(expi.trim)*sort(sdlevels[sdlevels!=0])
    expbins<-cut(expi,breaks=c(-Inf,cutoffs,Inf), labels=sort(sdlevels))
    names(expbins)=names(expi)
    binlevels=as.numeric(levels(expbins))
    for (bin in binlevels){
      if (bin>0){
        tmp.out=c(rownames(expMatrix)[i],paste("UP","_Level",abs(bin),sep=""),names(expbins[which(expbins %in% binlevels[binlevels>=bin])]))
      }else if(bin<0){
        tmp.out=c(rownames(expMatrix)[i],paste("DN","_Level",abs(bin),sep=""),names(expbins[which(expbins %in% binlevels[binlevels<=bin])]))
      }else{
        tmp.out=c()
      }
      if (length(tmp.out)>6){
        write(tmp.out,file=outfilepath,append=TRUE,sep="\t",ncolumns = length(tmp.out))  
      }
    }
  } 
}
#binarize expression features based on quartile transformed data
library(singscore)
binarize.expfeature.quartile<-function(expfile,levels=c(-6:-1,1:6),method=c("mean_sd","median_mad"),sdpercfilter=0.2,predefine.cutoff=NULL){#predefine.cutoff: file name for predefined cut offs for all genes
  #levels should not have 0, and should be symmetric, max number of levels=c(-10:-1,1:10), DO NOT EXCEED 10
  method <- match.arg(method)
  if (!grepl("\\.tsv",expfile)){
    print ("expfile must have extension of .tsv")
  }
  if (typeof(predefine.cutoff)!="NULL"){
    cutoffData<-read.delim(predefine.cutoff,row.names=1,stringsAsFactors = F,check.names = F,header = T,sep="\t")
    levels=colnames(cutoffData)
  }else{
    cutfile<-paste(sub("(.*?).tsv", "\\1", expfile),".",method,".",length(levels),"cutoffs.tsv",sep="")
    writeLines(paste(c("Gene",levels),collapse='\t'),cutfile)
  }
  gmtfile<-paste(sub("(.*?).tsv", "\\1", expfile),".",method,".",length(levels),"qrlevels.gmt",sep="")
  suppressWarnings(file.remove(gmtfile))
  expData<-read.table(expfile,stringsAsFactors = F,row.names = 1, check.names = F, header = T,sep = "\t",na.strings = c("", "NA"))
  if (length(grep("_1$", colnames(expData)))>0){
    expData<-expData[, -grep("_1$", colnames(expData))]#remove duplcate columns (Duplicated column names deduplicated: '1503362' => '1503362_1' [845])  
  }
  rownames(expData)<-gsub("\\.\\d+$", "", rownames(expData)) #remove .\d+ in ensemble_gene
  expMatrix<-rankGenes(expData)/nrow(expData)
  if (typeof(predefine.cutoff)=="NULL"){
    exp.sd <- apply(expMatrix, 1, sd)
    expMatrix<-expMatrix[which(exp.sd > quantile(exp.sd, sdpercfilter)), ]  
  }
  samples=colnames(expMatrix)
  writeLines(paste(c(paste("#Overlapping Expression Bins Generated Based On Cutoffs created using Trimmed (0.1) Mean + 1/2/3*SD based on log2 transformed fold change for quantile normalized TPM",sep=""),samples),collapse='\t'),gmtfile)
  for (i in 1:nrow(expMatrix)){
    expi=as.numeric(expMatrix[i,])
    names(expi)=names(expMatrix[i,])
    if (typeof(predefine.cutoff)!="NULL"){
      if (!rownames(expMatrix)[i] %in% rownames(cutoffData)){
        next
      }
      cutoffs=as.numeric(cutoffData[rownames(expMatrix)[i],])
    }else if (method=="mean_sd"){
      expi.trim=Trim(expi,trim=0.1)
      if (sd(expi.trim)==0){
        next
      }
      cutoffs=mean(expi.trim)+levels*sd(expi.trim)
      write(c(rownames(expMatrix)[i],cutoffs),file=cutfile,append=TRUE,sep="\t",ncolumns = length(cutoffs)+1)  
    }else if (method=="median_mad"){
      if (mad(expi)==0){
        next
      }
      cutoffs=median(expi)+levels*mad(expi)
      write(c(rownames(expMatrix)[i],cutoffs),file=cutfile,append=TRUE,sep="\t",ncolumns = length(cutoffs)+1)  
    }
    expbins<-cut(expi,breaks=c(-Inf,cutoffs,Inf), labels=c(levels[levels<0],0,levels[levels>0]))
    names(expbins)=names(expi)
    binlevels=as.numeric(levels(expbins))
    for (bin in binlevels){
      if (bin>0){
        tmp.out=c(rownames(expMatrix)[i],paste("UP_Level",bin,sep=""),names(expbins[which(expbins %in% binlevels[binlevels>=bin])]))
      }else if (bin<0){
        tmp.out=c(rownames(expMatrix)[i],paste("DN_Level",abs(bin),sep=""),names(expbins[which(expbins %in% binlevels[binlevels<=bin])]))
      }else{
        tmp.out=c()
      }
      if (length(tmp.out)>6){
        write(tmp.out,file=gmtfile,append=TRUE,sep="\t",ncolumns = length(tmp.out))  
      }
    }
  } 
}
#################################################################################   
###########Compute Feature Redundancy New Version v2021Sep Clustering based#####    
#################################################################################   
#fast calculate TCGA jaccard index matrix for genotypes of a subject and generate precalculation files
calSimilarity_subject_tcga<-function(subject.call,tcga.genotype.list,subject.genotype.list,method="ochiai",cutoff=0.1,by.cluster=TRUE,cutTree.method=c("quantile","dynamic"),dynamic.method="tree",minClusterSize,split.depth=2,visualize.cluster=FALSE,log.dir){ #method=="ochiai" or "jaccard"
  match.arg(cutTree.method)
  tmp.genoCell<-findConcepts(subject.call,subject.genotype.list)
  tmp.genoCell<-tmp.genoCell[tmp.genoCell %in% names(tcga.genotype.list)]
  tmp.tcga.genoList<-tcga.genotype.list[tmp.genoCell]
  if (length(tmp.tcga.genoList)==0){
    tmp.out<-c(0)
  }else{
    tmp.genolist.bin<-genolist2matrix(tmp.tcga.genoList)
    if (method=="ochiai"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "1-J/sqrt(A*B)"))
    }else if (method=="phi"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "((A-J)*(B-J)-J*(P-A-B+J))/sqrt(A*B*(P-A)*(P-B))"))
    }else if(method=="jaccard"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "(A+B-2*J)/(A+B-J)"))
    }else if(method=="raup-crick"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "1-phyper(J-1, A, P-A, B)"))
    }else if (method=="kulczynski"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "1-(J/2)*(1/A+1/B)"))
    }else{
      print ("method must be set as one of the following (pleae match case): ochiai,phi,jaccard,raup-crick, or kulczynski")
    }
    if (method=="phi"){
      sim.matrix=-dist.matrix  
    }else{
      sim.matrix=1-dist.matrix
    }
    sim.matrix <- ifelse(sim.matrix<cutoff,0,sim.matrix)
    if (by.cluster==TRUE & nrow(dist.matrix)>20){
      cluster.feature=fastcluster::hclust(as.dist(dist.matrix), method = "ward.D2")
      #groups=cutree(cluster.feature,h=Rfast::nth(unique(cluster.feature$height), 5, descending = T))
      if (cutTree.method=="quantile"){
        quantile.tree <- quantile(unique(cluster.feature$height), probs = seq(0, 1, 0.0001))
        groups=cutree(cluster.feature,h=quantile.tree["99.95%"]) 
      }else if (cutTree.method=="dynamic"){
        log <- capture.output({
          groups=cutreeDynamic(cluster.feature,minClusterSize = minClusterSize,method = dynamic.method,distM = dist.matrix,deepSplit = split.depth)
          names(groups)=cluster.feature$labels
        }) 
      }
      write(c(paste(length(unique(groups)),"feature groups detected for subject",subject.call)),file=paste0(log.dir,"/preCal.cluster.number.log"),append=TRUE)
      if (visualize.cluster==TRUE){
        dist.matrix.subject<-as.matrix(designdist(t(tmp.genolist.bin), method = "1-J/sqrt(A*B)"))
        cluster.subject=fastcluster::hclust(as.dist(dist.matrix.subject), method = "ward.D2")
        jpeg(paste0(log.dir,"/",gsub(" ",".",gsub(":","-",Sys.time())),".SUBJECT",subject.call,".cluster.heatmap.jpg"), width = 5000, height =5000)
        pheatmap(tmp.genolist.bin, cluster_rows=cluster.feature, cluster_cols=cluster.subject, cutree_rows=length(unique(groups)), show_rownames = FALSE, show_colnames = FALSE)
        graphics.off()
        jpeg(paste0(log.dir,"/",gsub(" ",".",gsub(":","-",Sys.time())),".SUBJECT",subject.call,".cluster.cormap.jpg"), width = 5000, height =5000)
        quantile.range <- quantile(sim.matrix, probs = seq(0, 1, 0.01))
        myBreaks <- seq(quantile.range["10%"], quantile.range["90%"], 0.01)
        myColor  <- colorRampPalette(c("skyblue", "white", "red"))(length(myBreaks) - 1)
        pheatmap(sim.matrix, cluster_rows=cluster.feature, cluster_cols=cluster.feature, breaks=myBreaks,color=myColor,cutree_rows=length(unique(groups)), cutree_cols=length(unique(groups)), show_rownames = FALSE, show_colnames = FALSE)
        graphics.off()
      }
      sum.similarity=c()
      for(i in 1:ncol(sim.matrix)){
        featurei=colnames(sim.matrix)[i]
        groupi=groups[featurei]
        groupfeatures=names(groups)[groups==groupi]
        sumi=sum(sim.matrix[groupfeatures,i])
        if (sumi==0){
          stop(paste("subject:",subject.call,"col:",i,"feature:",featurei,"group:",groupi,"sumi=0"))
        }
        sum.similarity=setNames(c(sum.similarity, sumi), c(names(sum.similarity), featurei))
      }
    }else{
      sum.similarity<-colSums(sim.matrix)
    }
    tmp.out=sum.similarity
  }
  return(tmp.out)
}
#precalculate dataset genotype redundacy based on TCGA and genotype lists
batch_calSimilarity_genolist<-function(subject.genotype.list,REF.genotype.list,subject.preCalfile,minTCGA=20,method="ochiai",cutoff=0.1,by.cluster=TRUE,cutTree.method=c("quantile","dynamic"),dynamic.method="tree",minClusterSize=40,split.depth=2,log.dir,visualize.cluter=FALSE){
  match.arg(cutTree.method)
  subject.id<-unique(unlist(subject.genotype.list))
  suppressWarnings(file.remove(subject.preCalfile))
  lostperc<-length(names(subject.genotype.list)[!names(subject.genotype.list) %in% names(REF.genotype.list)])/length(names(subject.genotype.list))
  header=paste0("#Precalculations of genotype redundancies: ","minTCGA=",minTCGA,", method=",method,", cutoff=",cutoff,", by.cluster=",by.cluster,", cutTree.method=",cutTree.method,", dynamic.method=",dynamic.method,", minClusterSize=",minClusterSize,", split.depth=",split.depth,". Loss of genotype from TCGA: ",100*lostperc,"%","\n")
  feature.redund=list()
  for (i in 1:length(subject.id)){
    if (i %% 100==0){
      print (paste("Processed",i,"subjects"))
      feature.redund[[subject.id[i]]]<-calSimilarity_subject_tcga(subject.call=subject.id[i],tcga.genotype.list=REF.genotype.list,subject.genotype.list=subject.genotype.list,method=method,cutoff=cutoff,by.cluster=by.cluster,cutTree.method=cutTree.method,dynamic.method=dynamic.method,minClusterSize=minClusterSize,split.depth=split.depth,visualize.cluster=visualize.cluter,log.dir = log.dir)
    }else{
      feature.redund[[subject.id[i]]]<-calSimilarity_subject_tcga(subject.call=subject.id[i],tcga.genotype.list=REF.genotype.list,subject.genotype.list=subject.genotype.list,method=method,cutoff=cutoff,by.cluster=by.cluster,cutTree.method=cutTree.method,dynamic.method=dynamic.method,minClusterSize=minClusterSize,split.depth=split.depth,visualize.cluster=F,log.dir = log.dir)
    }
  }
  result=list(feature.redundancy=feature.redund,feature.redundancy.parameters=header)
  save(result,file=subject.preCalfile)
}
##########################################################
#########   Feature Selection Modules          ###########
##########################################################
###calculate weight based on similarity or correlation statistics
calCorrelation<-function(target.list,compare.list,list.all,confound.factor=NULL,confound.cut=0.4){ #"pearson method"
  if (typeof(confound.factor)=="NULL"){
    myData <- data.frame(target.list=ifelse (list.all %in% target.list,1,0),compare.list=ifelse (list.all %in% compare.list,1,0))
    pearson=cor.test(myData$target.list,myData$compare.list,method="pearson")
    return(c(pearson$estimate,pearson$estimate,pearson$p.value,length(target.list),length(compare.list))) 
    }else{
    myData <- data.frame(target.list=ifelse (list.all %in% target.list,1,0),compare.list=ifelse (list.all %in% compare.list,1,0),confound.factor=confound.factor[match(list.all,names(confound.factor))])
    pearson=cor.test(myData$target.list,myData$compare.list,method="pearson")
    pearson.confound=cor.test(myData$compare.list,as.numeric(as.character(myData$confound.factor)),method="pearson")
    if (abs(pearson.confound$estimate)>confound.cut){
      return(c(0,0,1,length(target.list),length(compare.list)))
    }else{
      return(c(pearson$estimate,pearson$estimate,pearson$p.value,length(target.list),length(compare.list))) 
    }
  }
}
###calculate weight based on logistic regression statistics with consideration of confound factor
calLogReg<-function(target.list,compare.list,list.all,confound.factor=NULL,confound.pcut=0.01){
  if (typeof(confound.factor)=="NULL"){
    myData <- data.frame(target.list=ifelse (list.all %in% target.list,1,0),compare.list=ifelse (list.all %in% compare.list,1,0))
    LogReg=glm(formula = target.list ~ compare.list, family = "binomial", data = myData)
    LogReg.sum=summary(LogReg)$coefficients["compare.list",]
    #return(c(PseudoR2(LogReg, which = NULL),PseudoR2(LogReg, which = NULL),LogReg.sum["Pr(>|z|)"],length(target.list),length(compare.list))) 
    return(c(LogReg.sum["Estimate"],LogReg.sum["Estimate"],LogReg.sum["Pr(>|z|)"],length(target.list),length(compare.list))) 
  }else{
    myData <- data.frame(target.list=ifelse (list.all %in% target.list,1,0),compare.list=ifelse (list.all %in% compare.list,1,0),confound.factor=confound.factor[match(list.all,names(confound.factor))])
    #method 1 multivariate unconfound method
#    LogReg=glm(formula = target.list ~ compare.list + confound.factor, family = "binomial", data = myData)
#    LogReg.sum=summary(LogReg)$coefficients["compare.list",]
    #return(c(LogReg.sum["Estimate"],LogReg.sum["Estimate"],LogReg.sum["Pr(>|z|)"],length(target.list),length(compare.list))) 
    #method 2 univariate unconfound method
    LogReg=glm(formula = target.list ~ compare.list, family = "binomial", data = myData)
    LogReg.sum=summary(LogReg)$coefficients["compare.list",]
    LogReg.cf=glm(formula = target.list ~ confound.factor, family = "binomial", data = myData)
    LogReg.sum.cf=summary(LogReg.cf)$coefficients["confound.factor1",]
    if (LogReg.sum.cf["Pr(>|z|)"]<confound.pcut){
      return(c(0,0,1,length(target.list),length(compare.list)))
    }else{
      return(c(LogReg.sum["Estimate"],LogReg.sum["Estimate"],LogReg.sum["Pr(>|z|)"],length(target.list),length(compare.list)))
    }
  }
}
#Calculate weight based on ochiai Index
CalWeight.Similarity<-function(list1,list2,list.all,method){ #method=="ochiai" or "jaccard"  or phi
  J=length(intersect(list1,list2))
  A=length(list1)
  B=length(list2)
  P=length(list.all)
  J1=J-1;P1=P-1
  if (J1<0){
    J1=0
  }
  if (method=="ochiai"){
    tmp.weight=J/sqrt(A*B)
    tmp.weight.1=J1/sqrt(A*B)
  }else if (method=="jaccard"){
    tmp.weight=J/(A+B-J)
    tmp.weight.1=J1/(A+B-J1)
  }else if (method=="phi"){
    tmp.weight=(J*(P-A-B+J)-(A-J)*(B-J))/sqrt(A*B*(P-A)*(P-B))
    tmp.weight.1=(J1*(P1-A-B+J1)-(A-J1)*(B-J1))/sqrt(A*B*(P1-A)*(P1-B))
  }else{
    print ("method must be ochiai, jaccard, or phi")
  }
  return(c(tmp.weight,tmp.weight.1,J,A,B))
}

#batch calculate weight for a target list and a compendia of comparing lists
batch_CalWeight_Similarity<-function(target.list,compare.list,trainset=NULL,method=c("logistic.regression","pearson","phi","ochiai","jaccard"),minsize=10,confound.factor=NULL){ #method=="logistic.regression" or "ochiai" or "jaccard" or "pearson", or "Kappa"
  method <- match.arg(method)
  if (typeof(trainset)!="NULL"){
    target.list=target.list[which(target.list %in% trainset)]
    compare.list<-lapply(compare.list,filter,filter=trainset,inlist=TRUE)
  }
  compare.list<-compare.list[which(sapply(compare.list,function(x) length(x)>=minsize))]
  target.list=as.character(target.list)
  list.all=unique(unlist(compare.list))
  if (method=="pearson"){
    target.result<-data.frame(matrix(ncol = 6, nrow = 0))
    for (i in 1:length(compare.list)){
      target.result[nrow(target.result) + 1,]<-c(names(compare.list)[i],calCorrelation(target.list=target.list,compare.list=compare.list[[i]],list.all=list.all,confound.factor=confound.factor))
    }
    colnames(target.result)<-c("Feature.List","Weight","Weight.1","pValue","Target.Size","Compare.Size")
    } else if (method=="logistic.regression"){
    target.result<-data.frame(matrix(ncol = 6, nrow = 0))
    for (i in 1:length(compare.list)){
      target.result[nrow(target.result) + 1,]<-c(names(compare.list)[i],calLogReg(target.list=target.list,compare.list=compare.list[[i]],list.all=list.all,confound.factor=confound.factor))
    }
    colnames(target.result)<-c("Feature.List","Weight","Weight.1","pValue","Target.Size","Compare.Size")
    }
    else{
    target.result<-data.frame(matrix(ncol = 6, nrow = 0))
    for (i in 1:length(compare.list)){
      target.result[nrow(target.result) + 1,]<-c(names(compare.list)[i],CalWeight.Similarity(list1=target.list,list2=compare.list[[i]],list.all=list.all,method=method))
    }
    colnames(target.result)<-c("Feature.List","Weight","Weight.1","Intersect","Target.Size","Compare.Size")
  }
  return(target.result)
}


##########################################################
#########  Calculate GenSig modules            ###########
##########################################################
#remove equivocal features in the filtered feature weight table
#sig.feature.weight must contain Feature.List column
#feature.list.train must only contain the subject used in training sets
filter.feature<-function(feature.weight,feature.list.train,weight.cut,low.weight.cut,p.cut=NULL,low.p.cut,feature.col=1,highlevel.minsize=10,filter.criteria=c(1,2)){ #filter 1: Genotypes identified as low-level OE/UE but not high-level OE/UE, 2:Genotypes identified as both increased and decreased activity of the same gene
  if (typeof(p.cut)!="NULL"){
    feature.weight$selected=ifelse(as.numeric(feature.weight$Compare.Size)>=minsize&as.numeric(feature.weight$Weight)>weight.cut&as.numeric(feature.weight$pValue)<p.cut,1,0)
    feature.weight$selected.lowcut=ifelse(as.numeric(feature.weight$Compare.Size)>=minsize&as.numeric(feature.weight$Weight)>low.weight.cut&as.numeric(feature.weight$pValue)<low.p.cut,1,0)
  }else if (typeof(p.cut)=="NULL"){
    feature.weight$selected=ifelse(as.numeric(feature.weight$Compare.Size)>=minsize&as.numeric(feature.weight$Weight)>weight.cut,1,0)
    feature.weight$selected.lowcut=ifelse(as.numeric(feature.weight$Compare.Size)>=minsize&as.numeric(feature.weight$Weight)>low.weight.cut,1,0)
  }else {
    stop ("please provide p.cut in the function")
  }
  colnames(feature.weight)[feature.col]="Feature.List"
  feature.weight$ENSG=sub("\\:.*$","",feature.weight$Feature.List)
  feature.weight$activity=ifelse(grepl("UP_Level1",feature.weight$Feature.List),1,
                                 ifelse(grepl("UP\\_",feature.weight$Feature.List),2,
                                        ifelse(grepl("DN\\_Level1",feature.weight$Feature.List),-1,
                                               ifelse(grepl("DN\\_",feature.weight$Feature.List),-2,0))))
  up1only=unique(feature.weight$ENSG[feature.weight$activity == 1 & feature.weight$selected==1][which(!feature.weight$ENSG[feature.weight$activity == 1 & feature.weight$selected==1] %in% feature.weight$ENSG[feature.weight$activity == 2 & feature.weight$selected.lowcut==1])])
  dn1only=unique(feature.weight$ENSG[feature.weight$activity == -1 & feature.weight$selected==1][which(!feature.weight$ENSG[feature.weight$activity == -1 & feature.weight$selected==1] %in% feature.weight$ENSG[feature.weight$activity == -2 & feature.weight$selected.lowcut==1])])
  feature.length=stack(lengths(feature.list.train))
  feature.length=feature.length[str_count(feature.length$ind,":")==1,]
  feature.length=separate(data = feature.length, col = ind, into = c("ENSG", "ALT"), sep = ":")
  highlevel.up=feature.length[grepl("UP_Level([2-9]|10)",feature.length$ALT) & feature.length$values>=highlevel.minsize,]
  highlevel.dn=feature.length[grepl("DN_Level([2-9]|10)",feature.length$ALT) & feature.length$values>=highlevel.minsize,]
  reverse.weight.cut=low.weight.cut
  ENSG.contraversal=unique(c(intersect(feature.weight$ENSG[feature.weight$selected==1 & feature.weight$activity>0],feature.weight$ENSG[feature.weight$Weight>reverse.weight.cut & feature.weight$activity<0]),
                           intersect(feature.weight$ENSG[feature.weight$Weight>reverse.weight.cut & feature.weight$activity>0],feature.weight$ENSG[feature.weight$selected==1 & feature.weight$activity<0])))
  feature.weight.selected=feature.weight[feature.weight$selected==1,] #this step is critical for the following filtering determination. do not delete.
  feature.weight.selected$filter=ifelse(feature.weight.selected$activity==1 & (feature.weight.selected$ENSG %in% intersect(up1only,unique(highlevel.up$ENSG))),1,
                                        ifelse(feature.weight.selected$activity == -1 & (feature.weight.selected$ENSG %in% intersect(dn1only, unique(highlevel.dn$ENSG))),1,
                                               ifelse(feature.weight.selected$ENSG %in% ENSG.contraversal,2,0)))
  print (paste0(sum(feature.weight.selected$filter==1)," Genotypes identified as low-level OE/UE but not high-level OE/UE are found to be significant"))
  print (paste0(sum(feature.weight.selected$filter==2)," Genotypes identified as both increased and decreased activity of the same gene are found to be significant"))
  print (paste0(sum(feature.weight.selected$filter %in% filter.criteria)," features will be filered from ",nrow(feature.weight.selected)," significant features"))
  feature.weight.selected=feature.weight.selected[!feature.weight.selected$filter %in% filter.criteria,]
  return(feature.weight.selected[,!(colnames(feature.weight.selected) %in% c("ENSG","activity","filter","selected","selected.lowcut"))])
}

#calculate GenSig scores based on binary clinical feature
cal.GenSig.similarity<-function(target.list,feature.list,feature.select=NULL,trainset,preCalmatrix=NULL,TCGA.feature.list,weight.preCal=NULL,minsize=10,weight.cut,low.weight.cut=0.05,p.cut=NULL,low.p.cut=1,power=1,root=0.5,ECNpenalty=1,method=c("logistic.regression","pearson","phi","ochiai","jaccard"),redun.method="ochiai",redun.cut=0.1,rm.overfit=FALSE,rm.equivocal=TRUE,normalize.GenSig=TRUE,highlevel.minsize=10,confound.factor=confound.factor,gensig.model=NULL){ #method=="ochiai" or "jaccard" or "pearson", or "Kappa"
  method <- match.arg(method)
  feature.list<-feature.list[which(sapply(feature.list,function(x) length(x)>=minsize))]
  if (typeof(weight.preCal)=="NULL"){
    target.weight=batch_CalWeight_Similarity(target.list=target.list,compare.list=feature.list,trainset=trainset,method=method,confound.factor=confound.factor) 
  }else{
    target.weight=weight.preCal
    if (method=="pearson" | method=="logistic.regression"){
      colnames(target.weight)<-c("Feature.List","Weight","Weight.1","pValue","Target.Size","Compare.Size")  
    }else if (typeof(p.cut)=="NULL"){
      colnames(target.weight)<-c("Feature.List","Weight","Weight.1","Intersect","Target.Size","Compare.Size") 
    }else{
      stop("p.cut should be assigned NULL if non-correlation statistics is used")
    }
  }
  weight.result=target.weight
  if (typeof(gensig.model)=="NULL"){
    if (rm.equivocal){#remove genes with both positive activity (i.e. upregulation) and negative activitty (i.e. downregulation) in significant genotypes and genes with high level (2|3) over or under expression genotypes with more than ksMinSize samples, but only level 1 genotypes are significant
      target.weight=filter.feature(feature.weight=target.weight,feature.list.train=lapply(feature.list,function(x) x[x %in% trainset]),weight.cut=weight.cut,low.weight.cut=low.weight.cut,p.cut=p.cut,low.p.cut=low.p.cut,highlevel.minsize=highlevel.minsize)
    }else{
      if (typeof(p.cut)!="NULL"){
        target.weight=target.weight[as.numeric(target.weight$Compare.Size)>=minsize&as.numeric(target.weight$Weight)>weight.cut&as.numeric(target.weight$pValue)<p.cut,]
      }else if (typeof(p.cut)=="NULL"){
        target.weight=target.weight[as.numeric(target.weight$Compare.Size)>=minsize&as.numeric(target.weight$Weight)>weight.cut,]
      }else {
        stop ("please provide p.cut in the function")
      }
    }    
  }else{
    target.weight=target.weight[target.weight$Feature.List %in% gensig.model$Feature.select,]
  }
  print (paste(nrow(target.weight),"significant features identified after filtering"))
  if (typeof(feature.select)!="NULL"){
    print(paste0("number of features provided by user for filtering: ", length(feature.select)))
    target.weight=target.weight[target.weight$Feature.List %in% feature.select,]
    print(paste0("number of significant features after matching to user selected features: ", nrow(target.weight)))
  }
  if (typeof(preCalmatrix)=="NULL"){
    if (typeof(TCGA.feature.list)=="NULL"){
      stop ("TCGA.genotype.list must be provided if preCalmatrix is not provided") 
    }
    feature.select.list<-feature.list[which(names(feature.list) %in% target.weight$Feature.List)]
    if (length(feature.select.list)==0){
      stop (paste0("feature.select.list = 0 feature.list length=",length(feature.list)," feature.weight length=",nrow(target.weight)))
    }
  }
  result<-data.frame(matrix(ncol = 3, nrow = 0))
  subjects.id=unique(unlist(feature.list))
  for (i in 1:length(subjects.id)){
    subject.id=subjects.id[i]
    if (typeof(preCalmatrix)=="NULL"){
      tmp.line <- calSimilarity_subject_tcga(subject.call=subject.id,tcga.genotype.list=TCGA.feature.list,subject.genotype.list=feature.select.list,cutoff=redun.cut,method=redun.method)
      if(length(tmp.line)==2){
        result[nrow(result) + 1,]<-c(subject.id,0,ifelse(subject.id %in% target.list,1,0))
        next
      }
    }else{
      tmp.line<-unlist(strsplit(preCalmatrix[subject.id,],"\t"))
      if(length(tmp.line)==2){
        result[nrow(result) + 1,]<-c(subject.id,0,ifelse(subject.id %in% target.list,1,0))
        next
      }
    }
    tmp.episilon<-as.data.frame(do.call(rbind, strsplit(tmp.line[3:length(tmp.line)],"@")),stringsAsFactors=FALSE)
    colnames(tmp.episilon)<-c("Feature.List","Epsilon")
    #ECN=sum(1/(as.numeric(tmp.episilon[,"Epsilon"])^root))
    tmp.episilon[,"Epsilon"]=(as.numeric(tmp.episilon[,"Epsilon"]))^root
    ECN=nrow(tmp.episilon)/mean(as.numeric(tmp.episilon[,"Epsilon"]),trim=0.3)
    tmp.data<-merge(tmp.episilon,target.weight,by.x="Feature.List",by.y="Feature.List")
    tmp.data$Epsilon=as.character(tmp.data$Epsilon)
    tmp.data$Feature.List=as.character(tmp.data$Feature.List)
    if (rm.overfit==TRUE){
      if (subject.id %in% target.list){
        tmp.data$Weight4Cal=tmp.data$Weight.1
      }else{
        tmp.data$Weight4Cal=tmp.data$Weight
      }
    }else{
      tmp.data$Weight4Cal=tmp.data$Weight
    }
    tmp.data=cbind(tmp.data,normWeight=(as.numeric(tmp.data[,"Weight4Cal"])^power)/as.numeric(tmp.data[,"Epsilon"]))
    gensig.score<-sum(as.numeric(tmp.data[,"normWeight"]))/(ECN^ECNpenalty)
    result[nrow(result) + 1,]<-c(subject.id,gensig.score,ifelse(subject.id %in% target.list,1,0))
    if(i %% 100==0){
      print(paste("Processed ",i," sample",sep=""))
    }
  }
  colnames(result)<-c("SampleID","GenSig","Target.List")
  if (normalize.GenSig==TRUE){
    result$GenSig=normalize(as.numeric(result$GenSig))
  }
  return(list(Weight=weight.result,GenSig=result,Feature.select=target.weight$Feature.List))
}

#batch calculate uniGenSig scores based on preCalculated weights
#subject.sensitive must only contain subjects in the training set
batchCal.GenSig.similarity<-function(subject.sensitive=NULL,feature.select=NULL,weightfile=NULL,trainset=NULL,outprefix,subject.genotype.list,preCalmatrix=NULL,tcga.genotype.list=NULL,subject.gensigdir,method=c("logistic.regression","pearson","phi","ochiai","jaccard"),redun.method,redun.cut,minsize=10,sen.weightcut,res.weightcut,sen.pcut=NULL,res.pcut=NULL,power=2,root=0.5,ECNpenalty=1,rm.equivocal,highlevel.minsize=10,confound.factor=NULL,by.cluster=FALSE,validationset=FALSE){ #p.cut should be assigned NULL if similarity index is used instead of correlation statistics
  method <- match.arg(method)
  outfile<-paste(subject.gensigdir,"/",outprefix,"_GenSig.xls",sep="")
  if (typeof(preCalmatrix)!="NULL"){
    preCalmatrix.sel=matrix(ncol = 1, nrow = 0)
    for (i in 1:nrow(preCalmatrix)){
      tmp.line<-unlist(strsplit(preCalmatrix[i,],"\t"))
      subjectID=as.character(tmp.line[1])
      preCalmatrix.sel=rbind(preCalmatrix.sel,preCalmatrix[i,])
      rownames(preCalmatrix.sel)[nrow(preCalmatrix.sel)]=subjectID
    }
    preCalmatrix=preCalmatrix.sel
  }
  if (typeof(weightfile)=="NULL"){
    weightfile<-paste(subject.gensigdir,"/",outprefix,"_Weight.xls",sep="")
    subject.nontarget=trainset[which(!trainset %in% subject.sensitive)]
    print ("weight file not found, calculating weights and sensitive GenSig scores")
    if (by.cluster==TRUE){
      target.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".target"),target.list=subject.sensitive,feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,minsize=minsize,weight.cut=sen.weightcut,p.cut=sen.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,confound.factor=confound.factor)
    }else{
      target.result=cal.GenSig.similarity(target.list=subject.sensitive,feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,minsize=minsize,weight.cut=sen.weightcut,p.cut=sen.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,confound.factor=confound.factor)
    }
    print ("calculating weights and resistant GenSig scores")
    if (method %in% c("pearson","logistic.regression","phi")){
      target.weight=target.result[["Weight"]]
      target.weight[2:3]=sapply(target.weight[2:3],function(x){as.numeric(as.character(x))})
      nontarget.weight=target.result[["Weight"]]
      nontarget.weight[2:3]<-sapply(nontarget.weight[2:3],function(x){-as.numeric(as.character(x))})
      if (by.cluster==TRUE){
        nontarget.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".nontarget"),target.list=c(),feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,weight.preCal=nontarget.weight,minsize=10,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize)
      }else{
        nontarget.result=cal.GenSig.similarity(target.list=c(),feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,weight.preCal=nontarget.weight,minsize=10,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize)
      }
    }else{
      if (by.cluster==TRUE){
        nontarget.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".nontarget"),target.list=subject.nontarget,feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,minsize=minsize,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize)
      }else{
        nontarget.result=cal.GenSig.similarity(target.list=subject.nontarget,feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,minsize=minsize,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize)
      }
      target.weight=target.result[["Weight"]]
      nontarget.weight=nontarget.result[["Weight"]]
    }
    weight.preCal<-merge(target.weight,nontarget.weight,by.x="Feature.List",by.y="Feature.List",suffixes = c(".sensitive",".resistant"))
    write.table(weight.preCal,file=weightfile,row.names=FALSE,col.names=TRUE,sep="\t")
  }else if (typeof(subject.sensitive)=="NULL"){
    if (by.cluster==TRUE){
      if (validationset==TRUE){
        R.file=sub("_Weight.xls","_GenSig.rda",weightfile)
        if (file.exists(R.file)){
          gensig.model=mget(load(R.file, envir=(tmp<- new.env())), envir=tmp)$R.result
          print ("calculating sensitive GenSig scores")
          target.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".target"),target.list=c(),feature.list=subject.genotype.list,
                                                        feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                                        weight.preCal=gensig.model$target_result$Weight,minsize=10,weight.cut=sen.weightcut,p.cut=sen.pcut,power=power,root=root,ECNpenalty=ECNpenalty,
                                                        method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,
                                                        gensig.model=gensig.model$target_result
          )
          print ("calculating resistant GenSig scores")
          nontarget.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".nontarget"),target.list=c(),feature.list=subject.genotype.list,
                                                           feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                                           weight.preCal=gensig.model$nontarget_result$Weight,minsize=10,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,
                                                           method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,
                                                           gensig.model=gensig.model$nontarget_result
          )
        }else{
          stop(paste(R.file,"not found"))
        }
      }else{
        print ("Weight file found, calculating sensitive GenSig scores")
        weight.preCal=read.delim(weightfile,stringsAsFactors = F,check.names = F,header = T,sep="\t")
        target.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".target"),target.list=c(),feature.list=subject.genotype.list,
                                                      feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                                      weight.preCal=weight.preCal[,1:6],minsize=10,weight.cut=sen.weightcut,p.cut=sen.pcut,power=power,root=root,ECNpenalty=ECNpenalty,
                                                      method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,
                                                      gensig.model=NULL
        )
        print ("weight file found, calculating resistant GenSig scores")
        nontarget.result=cal.GenSig.similarity.bycluster(heatmap.outprefix=paste0(subject.gensigdir,"/",outprefix,".nontarget"),target.list=c(),feature.list=subject.genotype.list,
                                                         feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                                         weight.preCal=weight.preCal[,c(1,7:11)],minsize=10,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,
                                                         method=method,redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,
                                                         gensig.model=NULL
        )
      }
    }else{
      if (validationset==TRUE){
        R.file=sub("_Weight.xls","_GenSig.rda",weightfile)
        if (file.exists(R.file)){
          gensig.model=mget(load(R.file, envir=(tmp<- new.env())), envir=tmp)$R.result
        }else{
          stop(paste(R.file,"not found"))
        }
        print ("calculating sensitive GenSig scores")
        target.result=cal.GenSig.similarity(target.list=c(),feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                            weight.preCal=gensig.model$target_result$Weight, minsize=10,weight.cut=sen.weightcut,p.cut=sen.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,
                                            redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,gensig.model=gensig.model$target_result)
        print ("calculating resistant GenSig scores")
        nontarget.result=cal.GenSig.similarity(target.list=c(),feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                               weight.preCal=gensig.model$nontarget_result$Weight,minsize=10,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,
                                               redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,gensig.model=gensig.model$nontarget_result)
      }else{
        weight.preCal=read.delim(weightfile,stringsAsFactors = F,check.names = F,header = T,sep="\t")
        print ("calculating sensitive GenSig scores")
        target.result=cal.GenSig.similarity(target.list=c(),feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                            weight.preCal=weight.preCal[,c(1:6)], minsize=10,weight.cut=sen.weightcut,p.cut=sen.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,
                                            redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize)
        print ("calculating resistant GenSig scores")
        nontarget.result=cal.GenSig.similarity(target.list=c(),feature.list=subject.genotype.list,feature.select=feature.select,trainset=trainset,preCalmatrix=preCalmatrix,TCGA.feature.list=tcga.genotype.list,
                                               weight.preCal=weight.preCal[,c(1,7:11)],minsize=10,weight.cut=res.weightcut,p.cut=res.pcut,power=power,root=root,ECNpenalty=ECNpenalty,method=method,
                                               redun.method=redun.method,redun.cut=redun.cut,rm.overfit=TRUE,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize)
      }
    }
  }else{
    stop("either weigthfile or outprefix must be specified")
  }
  if (by.cluster==TRUE){
    colnames(target.result[["GenSig"]])[-1]=paste0(colnames(target.result[["GenSig"]])[-1],":target")
    colnames(nontarget.result[["GenSig"]])[-1]=paste0(colnames(nontarget.result[["GenSig"]])[-1],":nontarget")
    merge.score<-merge(target.result[["GenSig"]],nontarget.result[["GenSig"]],by.x="SampleID",by.y="SampleID")
  }else{
    merge.score<-merge(target.result[["GenSig"]][-3],nontarget.result[["GenSig"]][-3],by.x="SampleID",by.y="SampleID")
    names(merge.score)=c("SUBJECT_ID","GeneSig.sensitive","GeneSig.resistant")
    merge.score$Feature.sensitive=rep(length(target.result[["Feature.select"]]),times=nrow(merge.score))
    merge.score$Feature.resistant=rep(length(nontarget.result[["Feature.select"]]),times=nrow(merge.score))
  }
  write.table(merge.score,file=outfile,row.names=FALSE,col.names=TRUE,sep="\t")
  R.result=list(target_result=target.result,nontarget_result=nontarget.result)
  save(R.result,file=paste0(sub(".xls","",outfile),".rda"))
}


##########################################################
#########  Calculate dGenSig modules          ############
##########################################################
##find division line based on training set
find.dline<-function(gensig.result){
  gensig.result$slope<-gensig.result[,2]/gensig.result[,3]
  max.finite=max(gensig.result$slope[is.finite(gensig.result$slope)])
  gensig.result$slope<-ifelse(is.finite(gensig.result$slope),gensig.result$slope,(gensig.result[,2]+1)*max.finite)
  gensig.result<-gensig.result[order(gensig.result$slope,decreasing = TRUE),]
  n.sen<-c(0)
  n.res<-c(0)
  for(j in 1:nrow(gensig.result)){
    if(gensig.result$label[j]=="sen"){
      n.sen<-append(n.sen,n.sen[j]+1)
    }else{
      n.sen<-append(n.sen,n.sen[j])
    }
    if(gensig.result$label[j]=="res"){
      n.res<-append(n.res,n.res[j]+1)
    }else{
      n.res<-append(n.res,n.res[j])
    }
  }
  n.sen<-n.sen[2:length(n.sen)]
  n.res<-n.res[2:length(n.res)]
  gensig.result$perc.sen<-n.sen/max(n.sen)
  gensig.result$perc.res<-n.res/max(n.res)
  gensig.result$youden<-gensig.result$perc.sen+1-gensig.result$perc.res
  mySlope<-gensig.result[gensig.result$youden==max(gensig.result$youden),"slope"][1]
  return(mySlope)
}
### calculate dGenSig score based on linear method
cal.dGenSig.trial2trial.math<-function(gensig.file,dgensig.file,subject.train=NULL,phenoData){
  if(!"label" %in% colnames(phenoData)){
    stop("phenoData must contain a column called \"label\" indicating whether the subject is \"sen\" or \"res\" or \"mid\" or \"NA\"")
  }
  result.plot<-read.delim(gensig.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
  result.plot=result.plot[result.plot$SUBJECT_ID %in% phenoData$SUBJECT_ID,]
  result.plot[,2]<-(result.plot[,2]-min(result.plot[,2]))/(max(result.plot[,2])-min(result.plot[,2]))
  result.plot[,3]<-(result.plot[,3]-min(result.plot[,3]))/(max(result.plot[,3])-min(result.plot[,3]))
  #calculate dGenSig score, the line should be y-mySlope*x=0
  result.merge<-merge(result.plot,phenoData,by.x="SUBJECT_ID",by.y="SUBJECT_ID",all=TRUE)
  if (typeof(subject.train)=="NULL" | nrow(result.plot)-length(subject.train)==0){
    warning (paste("train set not found for",gensig.file,", using ALL subjects as trainset"))
    result.plot$TEST_SET =0
    mySlope<-find.dline(result.merge[!is.na(result.merge$label),])
  }else{
    result.plot$TEST_SET =ifelse(result.plot$SUBJECT_ID %in% subject.train,0,1) 
    result.train=result.merge[(!is.na(result.merge$label)) & (result.merge$SUBJECT_ID %in% subject.train),]
    mySlope<-find.dline(result.train)
  }
  result.plot$dGensig<-(result.plot$GeneSig.sensitive-mySlope*result.plot$GeneSig.resistant)/sqrt(1+mySlope^2)
  result.plot$mySlope<-rep(mySlope,nrow(result.plot))
  write.table(result.plot,file=dgensig.file, sep="\t",col.names=TRUE,row.names=FALSE)
}

###Batch calculate dGenSig score based on training clinical trial dataset
batchCal.dGenSig.trial2trial.train <- function(gensigdir,test.matrix=NULL,phenoData,method=c("math","ridge")){ #phenoData must contain a column called "label" indicating whether the subject is "sen" or "res" or "mid" or "NA"
  method=match.arg(method)
  if (typeof(test.matrix)!="NULL"){
    if (ncol(test.matrix)!=nrow(phenoData)){
      stop("test.matrix column number must match phenoData row number")
    }
  }
  files<-list.files(path = gensigdir, pattern = "_GenSig.xls$")
  for (file in files){
    if (typeof(test.matrix)=="NULL"){
      subject.train=NULL
    }else{
      if (grepl("Fold\\d+",file)){
        trainset<-test.matrix[str_match(file, "Fold([0-9a-z]+)")[,2],] #please adjust the pattern match based on the genSig file name patterns
      }else if (grepl("trainset\\_",file)){
        trainset<-test.matrix[str_match(file, "trainset\\_([0-9a-z]+)\\_")[,2],] #please adjust the pattern match based on the genSig file name patterns
      }else{
        stop(paste(file,"format wrong",sep=" "))
      }
      subject.train=colnames(trainset)[trainset==0]
    }
    prefix<-sub("_GenSig.xls", "", file)
    print (paste("analyzing",file,sep=" "))
    if (method=="math"){
      cal.dGenSig.trial2trial.math(
        gensig.file=paste(gensigdir,"/",file,sep=""),
        dgensig.file=paste(gensigdir,"/",sub("(.*?)_GenSig.xls", "\\1", file),"_",method,"_dGenSig.xls",sep=""),
        subject.train=subject.train,
        phenoData=phenoData
      )      
    }else if (method=="ridge"){
      cal.dGenSig.trial2trial.ml(
        gensig.file=paste(gensigdir,"/",file,sep=""),
        dgensig.file=paste(gensigdir,"/",sub("(.*?)_GenSig.xls", "\\1", file),"_",method,"_dGenSig.xls",sep=""),
        subject.train=subject.train,
        phenoData=phenoData
      )      
    }else{
      stop(paste0("method not found for ",method))
    }
  }
}

cal.dGenSig.trial2trial.math.validation<-function(gensig.file,dgensig.file,mySlope){
  result.plot<-read.delim(gensig.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
  result.plot[,2]<-(result.plot[,2]-min(result.plot[,2]))/(max(result.plot[,2])-min(result.plot[,2]))
  result.plot[,3]<-(result.plot[,3]-min(result.plot[,3]))/(max(result.plot[,3])-min(result.plot[,3]))
  #calculate dGenSig score, the line should be y-mySlope*x=0
  result.plot$dGensig<-(result.plot$GeneSig.sensitive-mySlope*result.plot$GeneSig.resistant)/sqrt(1+mySlope^2)
  result.plot$mySlope<-rep(mySlope,nrow(result.plot))
  write.table(result.plot,file=dgensig.file, sep="\t",col.names=TRUE,row.names=FALSE)
}

batchCal.dGenSig.trial2trial.validation <- function(train.gensigdir,validation.gensigdir,method=c("math","ridge")){
  method <- match.arg(method)
  files<-list.files(path = validation.gensigdir, pattern = "_GenSig.xls$")
  for (file in files){
    fold.call<-as.numeric(unlist(str_extract_all(file, "(?<=Fold)\\d+"))[1])
    prefix<-sub("_GenSig.xls", "", file)
    print (paste("analyzing",file,sep=" "))
    if (method=="math"){
      train.dgensig=read.delim(paste(train.gensigdir,"/",sub("(.*?)_GenSig.xls", "\\1", file),"_",method,"_dGenSig.xls",sep=""),stringsAsFactors = F,check.names = F,header = T,sep="\t") 
      cal.dGenSig.trial2trial.math.validation(
        gensig.file=paste(validation.gensigdir,"/",file,sep=""),
        dgensig.file=paste(validation.gensigdir,"/",sub("(.*?)_GenSig.xls", "\\1", file),"_",method,"_dGenSig.xls",sep=""),
        mySlope=train.dgensig$mySlope[1])
    }else if (method=="ridge"){
      dgensig.model=mget(load(paste(train.gensigdir,"/",sub("(.*?)_GenSig.xls", "\\1", file),"_",method,"_dGenSig.rda",sep=""), envir=(tmp<- new.env())), envir=tmp)
      cal.dGenSig.trial2trial.ml.validation(
        gensig.file=paste(validation.gensigdir,"/",file,sep=""),
        dgensig.file=paste(validation.gensigdir,"/",sub("(.*?)_GenSig.xls", "\\1", file),"_",method,"_dGenSig.xls",sep=""),
        dgensig.model=dgensig.model$dGenSig.model
      ) 
    }
  }
}

##########################################################
#########  Benchmark dGenSig modules          ############
##########################################################
#benchmark dGenSig scores based on binary response data and/or survival data
library(scales)
batchbenchmark.dGenSig.trial <- function(gensig.dir, phenoData,event.col=NULL,time.col=NULL,by.cluster,validationset=FALSE){ ##phenoData must contain a column called "label" indicating whether the subject is "sen" or "res" or "mid" or "NA", and must only contain the subjects from the treatment arm to be tested
  if(!"label" %in% colnames(phenoData)){
    stop("phenoData must contain a column called \"label\" indicating whether the subject is \"sen\" or \"res\" or \"mid\" or \"NA\"")
  }
  files<-list.files(path = gensig.dir, pattern = "_dGenSig.xls$",full.names=TRUE)
  benchmark.df<-data.frame(matrix(ncol = 10, nrow = 0))
  for (dGenSig.file in files){
    if (grepl("Fold\\d+",dGenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(dGenSig.file, "(?<=Fold)\\d+"))[1])
    }else if (grepl("trainset\\_",dGenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(dGenSig.file, "(?<=trainset\\_)\\d+"))[1]) 
    }else{
      stop(paste(dGenSig.file,"format wrong",sep=" "))
    }
    print (paste("analyzing file",basename(dGenSig.file),sep=" "))
    result.plot=read.delim(dGenSig.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
    if (sum(is.na(result.plot$dGensig))>0){
      warning (paste(paste(result.plot$SUBJECT_ID[is.na(result.plot$dGensig)],collapse="/"),"do not have dGenSig scores and are excluded from benchmarking"))
      result.plot=result.plot[!is.na(result.plot$dGensig),]  
    }
    if (validationset==TRUE){
      result.plot$TEST_SET=rep(1,nrow(result.plot)) 
    }
    basecol=ncol(result.plot)
    result.plot<-merge(result.plot,phenoData,by.x="SUBJECT_ID",by.y="SUBJECT_ID")
    result.plot$senlabel<-ifelse(is.na(result.plot$label),NA,ifelse(result.plot$label=="sen","sen","other"))
    colnames(result.plot)[which(colnames(result.plot)=="GeneSig.sensitive")]="GenSig:sen"
    colnames(result.plot)[which(colnames(result.plot)=="GeneSig.resistant")]="GenSig:res"
    if (sum(result.plot$TEST_SET)<5){
     warning(paste0("There are only ",sum(result.plot$TEST_SET)," subject in test set for fold",fold.call,", will use ALL subjects as test set"))
     result.plot$TEST_SET=1
    }
    test.plot=result.plot[result.plot$TEST_SET==1&(!is.na(result.plot$senlabel)),]
    tmp.roc<-roc(test.plot$senlabel,test.plot$dGensig,levels=c("other","sen"),quiet=FALSE,direction = "<")
    opt.cut<-coords(tmp.roc, x="best", input="threshold", best.method="youden",best.weights = c(2, 0.2),transpose = TRUE)
    pdf(file=paste0(gensig.dir,"/Fold",fold.call,"_iGenSig_","benchmark.pdf"),width=5, height=5)
    if ("GenSig:sen" %in% colnames(result.plot) & "GenSig:res" %in% colnames(result.plot) & "slope" %in% colnames(result.plot)){
      result.plot=result.plot[order(-result.plot$TEST_SET,result.plot$label,result.plot$dGensig,decreasing=TRUE),]
      if (is.na(fold.call)){
        plot.gensig=ggplot(result.plot,aes(x=`GenSig:res`,y=`GenSig:sen`)) + geom_point(shape=16,size=5,aes(color=label)) + geom_abline(intercept = 0, slope = as.numeric(unique(result.plot$mySlope)), colour="black", size=0.3, linetype="dotted", alpha=0.8) + theme(text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank()) + scale_color_manual(values = c("sen" = "red", "res" = "deepskyblue", "mid" = "grey"))
      }else{
        plot.gensig=ggplot(result.plot,aes(x=`GenSig:res`,y=`GenSig:sen`)) + geom_point(shape=21,stroke=0.2,size=5,aes(color=as.factor(TEST_SET),fill=label)) + geom_abline(intercept = 0, slope = as.numeric(unique(result.plot$mySlope)), colour="black", size=0.3, linetype="dotted", alpha=0.8) + theme(text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank()) + scale_color_manual(breaks = c(0, 1), values=c(rgb(255, 255, 255, alpha=0, names = NULL, maxColorValue = 255), "black")) + scale_fill_manual(values = c("sen" = "red", "res" = "deepskyblue", "mid" = "grey"))
      }
      print (plot.gensig)
    }
    graph.roc=ggroc(tmp.roc,color="red")+
      geom_point(aes(x=opt.cut[2], y=opt.cut[3]), colour="black", size=1)+
      geom_text(aes(x=opt.cut[2], y=opt.cut[3],label=paste("Cutoff:",round(opt.cut[1],digits=2),"(",percent(opt.cut[2]),",",percent(opt.cut[3]),")",sep="")),hjust=0, vjust=-0.5,colour="black",size=6)+
      theme(text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank())
    print(graph.roc)
    ann_colors = list(
      label = c(sen="gold", res="royalblue4",mid="gray95"),
      dGensig=colorRampPalette(c("white", "deepskyblue"))(n=500)
    )
    if (by.cluster==TRUE){
      feature.matrix=data.frame(row.names=result.plot$SUBJECT_ID,result.plot[,grepl("GenSig\\d+",colnames(result.plot)),drop=FALSE])
    }else{
      feature.matrix=data.frame(row.names=result.plot$SUBJECT_ID,result.plot[,grepl("GenSig:",colnames(result.plot))]) 
    }
    feature.matrix=t(scale(feature.matrix,center=TRUE,scale=TRUE))
    quantile.range <- quantile(feature.matrix[feature.matrix!=0], probs = seq(0, 1, 0.01))
    paletteLength <- 50
    myColor <- colorRampPalette(c("skyblue", "white", "red"))(paletteLength)
    myBreaks <- c(seq(quantile.range["1%"], 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(quantile.range["99%"]/paletteLength, quantile.range["99%"], length.out=floor(paletteLength/2)))
    ann_subject=data.frame(row.names=result.plot$SUBJECT_ID,result.plot[,c("TEST_SET","label","dGensig")])
    ann_subject=ann_subject[order(ann_subject$TEST_SET,ann_subject$dGensig),]
    feature.matrix=feature.matrix[,rownames(ann_subject),drop=FALSE]
    #if (nrow(feature.matrix)>1){
      heatmap=pheatmap(feature.matrix, cluster_rows=TRUE, cluster_cols=FALSE, treeheight_row = 0, color=myColor, breaks=myBreaks,show_rownames = T, show_colnames = F,annotation_col = ann_subject,annotation_colors = ann_colors,legend=F)
      print (heatmap)      
    #}
    graphics.off()
    if (is.null(event.col)==FALSE){
      surv<-survival.analysis.trial(dGenSig=result.plot[result.plot$TEST_SET==1,],dGenSig.file=dGenSig.file,event.col=event.col+basecol-1,time.col=time.col+basecol-1,rho=0)
    }else{
      surv=c(NA,NA,NA,NA)
    }
    write.table(result.plot,file=paste(sub(".xls","",dGenSig.file),"_phenoData.xls",sep=""), sep="\t",col.names=TRUE,row.names=FALSE)
    if (sum(grepl("Feature.sensitive|Feature.resistant",colnames(test.plot)))<2){
    Feature.count=c(NA,NA)
    }else{
    Feature.count=c(test.plot$Feature.sensitive[1],test.plot$Feature.resistant[1])  
    }
    for(predictor in colnames(test.plot)[grepl("gensig",colnames(test.plot),ignore.case = T)]){
      if (grepl("res|nontarget",predictor)){
        tmp.roc<-roc(test.plot[,"senlabel"],test.plot[,predictor],levels=c("other","sen"),quiet=FALSE,direction = ">")
      }else{
        tmp.roc<-roc(test.plot[,"senlabel"],test.plot[,predictor],levels=c("other","sen"),quiet=FALSE,direction = "<")
      }
      benchmark.df[nrow(benchmark.df) + 1,]<-c(fold.call,predictor,length(result.plot$SUBJECT_ID),tmp.roc$auc,surv,Feature.count)
    }
    if ("GenSig:sen" %in% colnames(result.plot) & "GenSig:res" %in% colnames(result.plot)){
    logR.matrix=data.frame(row.names=test.plot$SUBJECT_ID,test.plot[,c("GenSig:sen","GenSig:res")],Class=ifelse(test.plot$label=="sen",1,0))
    model.train <- glm(Class ~ ., data = logR.matrix, family = "binomial", maxit = 100)
    logR.predictor=predict(model.train, newdata = logR.matrix[,!colnames(logR.matrix) %in% c("Class")], type = "response")
    logR.roc=roc(logR.matrix$Class,as.numeric(logR.predictor),levels=c("0","1"),quiet=FALSE,direction = "<")
    benchmark.df[nrow(benchmark.df) + 1,]<-c(fold.call,"logR.igensig",length(result.plot$SUBJECT_ID),logR.roc$auc,surv,Feature.count)
    }
  }
  colnames(benchmark.df)<-c("Fold","predictor","SUBJECT(n=)","AUC","Survival.optCut","HazardRatio","PValue","Spearman.Corr.CoxResidue","Feature.sensitive","Feature.resistant")
  if ("GenSig:sen" %in% benchmark.df$predictor){
    for(predictor in c("GenSig:sen","GenSig:res","dGensig","logR.igensig")){
      benchmark.df.tmp=benchmark.df[grepl(paste0("^",predictor,"$"),benchmark.df$predictor) & benchmark.df$Fold!=0,-c(1:2)] #fold=0 indicate no permutation
      benchmark.df=rbind(benchmark.df, data.frame(Fold="average",predictor=predictor,t(apply(benchmark.df.tmp, 2, function (x) mean(as.numeric(x))) ),check.names = F)) #colMedians(benchmark.df[,-1])
    } 
  }else{
    predictor="dGensig"
    benchmark.df.tmp=benchmark.df[grepl(paste0("^",predictor,"$"),benchmark.df$predictor),-c(1:2)]
    benchmark.df=rbind(benchmark.df, data.frame(Fold="average",predictor=predictor,t(apply(benchmark.df.tmp, 2, function (x) mean(as.numeric(x))) ),check.names = F)) #colMedians(benchmark.df[,-1])
  }
  write.table(benchmark.df,file=paste(gensig.dir,"/",Sys.Date(),".dGenSig_",sample(1:99, 1),"_benchmark.result.xls",sep=""),sep="\t", row.names=FALSE, col.names = TRUE)
}
#benchmark dGenSig scores based on binary response data via combining all folds
batchbenchmark.dGenSig.combinefolds.trial <- function(gensig.dir, phenoData,event.col=NULL,time.col=NULL){ ##phenoData must contain a column called "label" indicating whether the subject is "sen" or "res" or "mid" or "NA", and must only contain the subjects from the treatment arm to be tested
  if(!"label" %in% colnames(phenoData)){
    stop("phenoData must contain a column called \"label\" indicating whether the subject is \"sen\" or \"res\" or \"mid\" or \"NA\"")
  }
  files<-list.files(path = gensig.dir, pattern = "_dGenSig.xls$",full.names=TRUE)
  combine.plot<-data.frame(matrix(ncol = 6, nrow = 0))
  colnames(combine.plot)=c("SUBJECT_ID","GeneSig.sensitive","GeneSig.resistant","TEST_SET","dGensig","mySlope")
  for (dGenSig.file in files){
    if (grepl("Fold\\d+",dGenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(dGenSig.file, "(?<=Fold)\\d+"))[1])
    }else if (grepl("trainset\\_",dGenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(dGenSig.file, "(?<=trainset\\_)\\d+"))[1]) 
    }else{
      stop(paste(dGenSig.file,"format wrong",sep=" "))
    }
    result.plot=read.delim(dGenSig.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
    result.plot$dGensig=(result.plot$dGensig-mean(Trim(result.plot$dGensig,trim=0.1)))/sd(Trim(result.plot$dGensig,trim=0.1))
    result.plot$dGensig=normalize(result.plot$dGensig)
    combine.plot=rbind(combine.plot,result.plot[result.plot$TEST_SET==1,])
  }
  combine.plot<-merge(combine.plot,phenoData,by.x="SUBJECT_ID",by.y="SUBJECT_ID")
  combine.plot$senlabel<-ifelse(is.na(combine.plot$label),NA,ifelse(combine.plot$label=="sen","sen","other"))
  write.table(combine.plot,file = paste(gensig.dir,"/",basename(gensig.dir),".dGenSig",".combinetestsets.xls",sep=""), quote = F, row.names = F, col.names = T, sep="\t")
  combine.test=combine.plot[combine.plot$TEST_SET==1&(!is.na(combine.plot$senlabel)),]
  tmp.roc<-roc(combine.test$senlabel,combine.test$dGensig,levels=c("other","sen"),quiet=FALSE,direction = "<")
  print (tmp.roc$auc)
}

#benchmark GenSig scores based on binary response data and/or survival data
batchbenchmark.GenSig.trial <- function(gensig.dir, phenoData,event.col=NULL,test.matrix=NULL,time.col=NULL,validationset=FALSE){ ##phenoData must contain a column called "label" indicating whether the subject is "sen" or "res" or "mid" or "NA", and must only contain the subjects from the treatment arm to be tested
  if(!"label" %in% colnames(phenoData)){
    stop("phenoData must contain a column called \"label\" indicating whether the subject is \"sen\" or \"res\" or \"mid\" or \"NA\"")
  }
  files<-list.files(path = gensig.dir, pattern = "_GenSig.xls$",full.names=TRUE)
  benchmark.df<-data.frame(matrix(ncol = 7, nrow = 0))
  for (GenSig.file in files){
    if (grepl("Fold\\d+",GenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(GenSig.file, "(?<=Fold)\\d+"))[1])
      if (validationset==F){
        trainset<-test.matrix[str_match(GenSig.file, "Fold([0-9a-z]+)")[,2],] 
        subject.train=colnames(trainset)[trainset==0]
      }
    }else{
      stop(paste(GenSig.file,"format wrong",sep=" "))
    }
    print (paste("analyzing file",basename(GenSig.file),sep=" "))
    result.plot=read.delim(GenSig.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
    result.plot$dGensig=result.plot$GeneSig.sensitive
    if (validationset==F){
      result.plot$TEST_SET =ifelse(result.plot$SUBJECT_ID %in% subject.train,0,1) 
      if (sum(result.plot$TEST_SET==1)==0){
        warning ("not test subjects found, using ALL subject as testing set")
        result.plot$TEST_SET =1  
      }
    }else{
      result.plot$TEST_SET =1
    }
    basecol=ncol(result.plot)
    result.plot<-merge(result.plot,phenoData,by.x="SUBJECT_ID",by.y="SUBJECT_ID")
    result.plot$senlabel<-ifelse(is.na(result.plot$label),NA,ifelse(result.plot$label=="sen","sen","other"))
    test.plot=result.plot[result.plot$TEST_SET==1&(!is.na(result.plot$senlabel)),]
    tmp.roc<-roc(test.plot$senlabel,test.plot$dGensig,levels=c("other","sen"),quiet=FALSE,direction = "<")
    if (is.null(event.col)==FALSE){
      surv<-survival.analysis.trial(dGenSig=result.plot[result.plot$TEST_SET==1,],dGenSig.file=GenSig.file,event.col=event.col+basecol-1,time.col=time.col+basecol-1,rho=0)
    }else{
      surv=c(NA,NA,NA,NA)
    }
    write.table(result.plot,file=paste(sub(".xls","",GenSig.file),"_phenoData.xls",sep=""), sep="\t",col.names=TRUE,row.names=FALSE)
    benchmark.df[nrow(benchmark.df) + 1,]<-c(fold.call,length(result.plot$SUBJECT_ID),tmp.roc$auc,surv)
  }
  colnames(benchmark.df)<-c("Fold","SUBJECT(n=)","AUC","Survival.optCut","HazardRatio","PValue","Spearman.Corr.CoxResidue")
  benchmark.df=rbind(benchmark.df, data.frame(Fold="Average",t(colMeans(benchmark.df[,-1])),check.names = F))
  write.table(benchmark.df,file=paste(gensig.dir,"/",basename(gensig.dir),".GenSig_","benchmark.result.xls",sep=""),sep="\t", row.names=FALSE, col.names = TRUE)
}
#benchmark dGenSig scores based on binary response data and/or survival data
batchbenchmark.dGenSig.trial.stratified <- function(gensig.dir, phenoData,pivot,event.col=NULL,time.col=NULL){ ##phenoData must contain a column called "label" indicating whether the subject is "sen" or "res" or "mid" or "NA", and must only contain the subjects from the treatment arm to be tested
  if(!"label" %in% colnames(phenoData)){
    stop("phenoData must contain a column called \"label\" indicating whether the subject is \"sen\" or \"res\" or \"mid\" or \"NA\"")
  }
  files<-list.files(path = gensig.dir, pattern = "_dGenSig.xls$",full.names=TRUE)
  benchmark.df<-data.frame(matrix(ncol = 6, nrow = 0))
  for (dGenSig.file in files){
    if (grepl("Fold\\d+",dGenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(dGenSig.file, "(?<=Fold)\\d+"))[1])
    }else if (grepl("trainset\\_",dGenSig.file)){
      fold.call<-as.numeric(unlist(str_extract_all(dGenSig.file, "(?<=trainset\\_)\\d+"))[1]) 
    }else{
      stop(paste(dGenSig.file,"format wrong",sep=" "))
    }
    print (paste("analyzing file",basename(dGenSig.file),sep=" "))
    result.plot=read.delim(dGenSig.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
    basecol=ncol(result.plot)
    result.plot<-merge(result.plot,phenoData,by.x="SUBJECT_ID",by.y="SUBJECT_ID")
    result.plot$senlabel<-ifelse(is.na(result.plot$label),NA,ifelse(result.plot$label=="sen","sen","other"))
    test.plot=result.plot[result.plot$TEST_SET==1&(!is.na(result.plot$senlabel)),]
    tmp.roc<-roc(test.plot$senlabel,test.plot$dGensig,levels=c("other","sen"),quiet=FALSE,direction = "<")
    if (is.null(event.col)==FALSE){
      pdf(file=paste(sub("_dGenSig.xls","",dGenSig.file),"_dGenSig.sensitive.pdf",sep=""),width=4, height=10)
      result.plot$stratification=ifelse(result.plot$GeneSig.sensitive>median(result.plot$GeneSig.sensitive),1,0)
      survival.analysis.stratified(survData=result.plot[result.plot$TEST_SET==1,],subject.col=1,genotype.col=which(colnames(result.plot)=="stratification"),event.col=event.col+basecol-1,time.col=time.col+basecol-1,pivot.col=which(colnames(result.plot)==pivot),rho=0,predict.sensitive=FALSE)
      dev.off()
      pdf(file=paste(sub("_dGenSig.xls","",dGenSig.file),"_dGenSig.resistant.pdf",sep=""),width=4, height=10)
      result.plot$stratification=ifelse(result.plot$GeneSig.resistant>median(result.plot$GeneSig.resistant),1,0)
      survival.analysis.stratified(survData=result.plot[result.plot$TEST_SET==1,],subject.col=1,genotype.col=which(colnames(result.plot)=="stratification"),event.col=event.col+basecol-1,time.col=time.col+basecol-1,pivot.col=which(colnames(result.plot)==pivot),rho=0,predict.sensitive=FALSE)
      dev.off()
      surv=c(NA,NA,NA)
    }else{
      surv=c(NA,NA,NA)
    }
    write.table(result.plot,file=paste(sub(".xls","",dGenSig.file),"_phenoData.xls",sep=""), sep="\t",col.names=TRUE,row.names=FALSE)
    benchmark.df[nrow(benchmark.df) + 1,]<-c(fold.call,length(result.plot$SUBJECT_ID),tmp.roc$auc,surv)
  }
  colnames(benchmark.df)<-c("Fold","SUBJECT(n=)","AUC","Survival.optCut","HazardRatio","PValue")
  write.table(benchmark.df,file=paste(gensig.dir,"/benchmark.result.xls",sep=""),sep="\t", row.names=FALSE, col.names = TRUE)
}
#########  Survival analysis module for dGenSig
library(maxstat)
survival.analysis.trial <- function(dGenSig,dGenSig.file,event.col=8,time.col=9,rho=0){
  colnames(dGenSig)[c(event.col,time.col)]=c("event","time")
  dGenSig=dGenSig[!is.na(dGenSig$event),]
  dGenSig$time=as.numeric(dGenSig$time)
  if (sum(dGenSig$dGensig!=0)>=3){
    condMC=maxstat.test(Surv(time, event) ~ dGensig, data = dGenSig, smethod = "LogRank", pmethod="condMC", B = 10000)
    cutoff=condMC$estimate
  }else{
    print (paste0(sum(dGenSig$dGensig !=0)," subjects have nonzero dGenSig scores, will use 0 as cutoff for", pdf.file))
    cutoff=0
  }
  methGroups = ifelse(dGenSig$dGensig>cutoff,"1High","2Low")
  hnSurv <- survfit(Surv(time,event)~as.factor(methGroups), 
                    data = dGenSig, conf.int=.98, se.fit=T, type=c("kaplan-meier"),
                    error=c("g"),conf.type=c("plain"),conf.lower=c("usual"))
  par(mfrow=c(2,1))
  # plot the functional forms
  pdf(file=paste(sub(".xls","",dGenSig.file),"_survival.pdf",sep=""),width=4, height=4)
  fittx<-coxph(Surv(time,event)~dGensig,data=dGenSig)
  targetval<-dGenSig[names(resid(fittx)),"dGensig"]
  plot(targetval,
       resid(fittx),cex=0.5, xlab=paste("Cutoff=",round(cutoff,digits=3)),ylab='Excess relapse or death'
  )
  cor.tmp=cor.test(targetval,as.numeric(resid(fittx)),method="pearson")
  smooth<-lowess(targetval,resid(fittx),iter=0)
  lines(smooth)
  abline(v=cutoff,lty=2)
  hr<-1/exp(coxph(Surv(time,event)~as.factor(methGroups), data = dGenSig)$coef)
  sdf <- survdiff(Surv(time,event)~as.factor(methGroups),data = dGenSig, rho=rho)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  plot(hnSurv,col=c("red","blue"), bty="l", lwd=2, las=1, xlab=paste("Coxph HR=",round(hr,digits=2),"Survdiff P=",round(p.val,digits=3), "rho=",rho))
  print (paste(sub(".xls","",dGenSig.file),"_survival.pdf"," created",sep=""))
  graphics.off()
  return(c(cutoff,hr,p.val,cor.tmp$estimate))
}

#find the signature pathways in the selected significant features provided in weight.tables
#weight.tables object must be a list of weight tables filtered based in significant criteria with first column named as "Feature.List"
#out.prefix: the path and prefix of output files
find.signature.pathways<-function(weight.tables,out.prefix,rm.equivocal=FALSE,feature.list.train,highlevel.minsize=20,concept.gmt,concept.preCal,pathway.list,CSEA.module,biomart.file=NULL){
  source(CSEA.module)
  if (typeof(biomart.file)=="NULL"){
    library('biomaRt')
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) 
  }else{
    mart=read.delim(biomart.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
  }
  concept.list=read_concepts(concept.gmt,min=10)
  concept.preCalmatrix<- as.matrix(unlist(readLines(concept.preCal)[-1]),col=1)
  for (j in 1:length(weight.tables)){
    weight.table.select=weight.tables[[j]]
    weight.table.type=names(weight.tables)[j]
    weight.table.select$ENSG=sub("\\:.*$","",weight.table.select$Feature.List)
    if (rm.equivocal==TRUE){
      weight.table.select=rm.equivocal.feature(sig.feature.weight=weight.table.select,feature.list.train=feature.list.train,highlevel.minsize=highlevel.minsize)
    }
    activMutations=weight.table.select$ENSG[grepl("Hotspot_Mutations|Mutation_Hotspots",weight.table.select$Feature.List)]
    weight.table.select$activity=ifelse(grepl("(Hotspot_Mutations|Mutation_Hotspots|UP\\_)",weight.table.select$Feature.List),1,
                                        ifelse(grepl("DN\\_",weight.table.select$Feature.List),0,
                                               ifelse(grepl("Nonsynonymous_Mutations",weight.table.select$Feature.List)&(weight.table.select$ENSG %in% activMutations),1,
                                                      ifelse(grepl("Nonsynonymous_Mutations",weight.table.select$Feature.List),0,NA)
                                               )))
    weight.table.select$equivocal=ifelse(weight.table.select$ENSG %in% intersect(unique(weight.table.select$ENSG[weight.table.select$activity==0]),unique(weight.table.select$ENSG[weight.table.select$activity==1])),1,0)
    print (paste("ratio of equivocal genes:",length(unique(weight.table.select$ENSG[weight.table.select$equivocal==1]))/length(unique(weight.table.select$ENSG)),sep=""))
    for (i in -1:1){
        if (i==-1){
        ENSG.selected=unique(weight.table.select$ENSG)
        }else{
        ENSG.selected=unique(weight.table.select$ENSG[weight.table.select$activity==i]) 
        }
        if (typeof(biomart.file)=="NULL"){
          target.list<-getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol"),values=ENSG.selected,mart= mart)$hgnc_symbol
        }else{
          target.list=mart$`Approved symbol`[match(ENSG.selected,mart$`Ensembl gene ID`)]
          target.list=unique(target.list[!is.na(target.list)])
        }
        uniConSig=cal.uniConSig(target.list=target.list,feature.list=concept.list,preCalmatrix=concept.preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai",rm.overfit=FALSE)
        if(typeof(uniConSig)=="NULL"){
           print (paste0("The target list is functionally heterogenous for ",weight.table.type," , CSEA pathway analysis cannot be carried out"))
           next
         }else if(nrow(uniConSig)<2000){
           print (paste0("number of genes having uniConSig scores are ",nrow(uniConSig),", cannot perform CSEA for",ifelse(i==-1,"UPDN",ifelse(i==1,"UP","DN")), " genes for predicting ",weight.table.type))
           next
         }
        CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list=pathway.list,p.cut=0.05,minsize=5,min.numOnList=5,transformNegWeight=FALSE)
        if(!exists("CSEA.result")){
          print (paste0("No significant ", ifelse(i==-1,"UPDN",ifelse(i==1,"UP","DN")), " pathways detected by CSEA in the genomic signature for predicting ",weight.table.type))
          next
        }
        disambiguate<-disambiguation.CSEA(GEA.result=CSEA.result,uniConSig.result=uniConSig,compare.list=pathway.list,upPathways=NULL,topn=min(c(30,nrow(CSEA.result))),p.cut=0.01)
        assoc <- pathwayAssociation(topPathway=disambiguate[[1]]$Compare.List[1:min(c(30,nrow(disambiguate[[1]])))],compare.list=pathway.list,feature.list=concept.list,preCalmatrix=concept.preCalmatrix,minsize=10,rm.overfit=FALSE)
        pathway.result=list(CSEA=CSEA.result,disambiguate=disambiguate,associations=assoc)
        save(pathway.result, file = paste(out.prefix,"_",weight.table.type,"_",ifelse(i==-1,"UPDN",ifelse(i==1,"UP","DN")),"_Pathway.RData",sep=""))
        if (nrow(pathway.result$associations)>1 & sum(pathway.result$associations)>=sum(diag(pathway.result$associations))){
          pdfsize=max(10,ncol(pathway.result$associations))
          pdf(file=paste(out.prefix,"_",weight.table.type,"_",ifelse(i==-1,"UPDN",ifelse(i==1,"UP","DN")),"_Pathway.pdf",sep=""),width=pdfsize, height=pdfsize)
          pathway.plot=pathway.corplot(pathway.result$associations)
          graphics.off()
        }else{
          print (paste0("only one significant pathway or there is no association between pathways detected in the diambiguated result"))
        }
      }
  }
}

find.pathway.feature<-function(select.pathways,HGNC.filter=NULL,concept.gmt,concept.preCal,pathway.gmt,CSEA.module,minSize=30){
  source(CSEA.module)
  concept.list=read_concepts(concept.gmt,min=10)
  concept.preCalmatrix<- as.matrix(unlist(readLines(concept.preCal)[-1]),col=1)
  pathway.list=c(read_concepts(pathway.gmt,min=10))
  select.pathway.list=pathway.list[select.pathways]
  if (length(select.pathway.list)==0){
    stop("selected pathways not found in the pathway gene set provided, please provide the pathway gene set containing the pathway names selected")
  }
  select.pathway.genes=unique(unlist(select.pathway.list))
  if (length(select.pathway.genes)<minSize){
    print (paste0("only ",length(select.pathway.genes)," genes found in the selected pathways. A minimal of ",minSize," is required"))
    return(list())
  }else{
    if (typeof(HGNC.filter)!="NULL"){
      select.pathway.genes=select.pathway.genes[select.pathway.genes %in% HGNC.filter]
    }
    select.uniConSig=cal.uniConSig(target.list=select.pathway.genes,feature.list=concept.list,concept.preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="ochiai",rm.overfit=FALSE)
    test.roc<-roc(select.uniConSig$Target.List,select.uniConSig$uniConSig,levels=c(0,1),quiet=FALSE,direction = "<")
    opt.cut<-coords(test.roc, x="best", input="threshold", best.method="youden",best.weights = c(1, 0.5),transpose=TRUE)
    pathway.related.genes=select.uniConSig$subjectID[select.uniConSig$uniConSig>opt.cut["threshold"] ]
    result.list=list(compute.pathway.genes=pathway.related.genes,selected.pathway.genes=select.pathway.genes,select.pathways=select.pathways,pathway.uniConSig=select.uniConSig,uniConSig.cutoff=opt.cut)
    return(result.list)  
  }
}

