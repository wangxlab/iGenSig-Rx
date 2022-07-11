# How to use iGenSig-Rx R modules to build multi-omics models for predicting therapeutic responses based on GDSC dataset

## iGenSig-Rx: an integral genomic signature based white-box method for modeling clinical therapeutic responses using genome-wide sequencing data

•	The current version of iGenSig-Rx is beta 3.2.3 (July 11th, 2022). 

•	The iGenSig-Rx was built on R version 4.1.0 and tested on Linux and Windows environment. 

## A. Introduction

•	Here we developed an integral genomic signature (iGenSig-Rx) analysis as a new class of transparent, interpretable, and resilient modeling methods for big data-based precision medicine with improved cross-dataset applicability and tolerance to sequencing bias. 

•	We define genomic features that significantly correlate with a clinical phenotype (such as therapeutic response) as genomic correlates, and an integral genomic signature as the integral set of redundant genomic correlates for a given clinical phenotype such as therapeutic response.

•	We postulate that the redundant genomic features, which are usually eliminated through dimensionality reduction or feature removal  during multi-omics modeling, may help overcome the sequencing bias. 

•	Using genomic dataset for chemical perturbations, we developed the iGenSig-Rx models predicting clinical response to HER2-targeted therapy and tested the cross-dataset performance of selected models based on independent multi-omics datasets for breast cancer HER2-positive patients assessing their therapeutic responses.

## B. How to download R code and data files to run iGenSig-Rx

• Download all files through git to your local directory. 

$ git clone https://github.com/wangxlab/iGenSig-Rx

If you have downloaded all necessary data files, you can find these 4 directories; 

1	“./GenotypePrecalmatrixData

You can find rds files for genotype and precalculated feature redundancy data of 4 datasets (CALGB, ACOSOG, NOAH, and TCGA).

2	“./PhenotypeData” directory

You can find Phenotype data files (.tsv file) of 3 datasets (CALGB, ACOSOG, and NOAH).

3.	“./TestsetAnnotationData” directory

You can find annotation files containing the subject IDs of permutated CALGB test sets (.tsv file). 


## C. Installing R packages

•	Find "iGenSigRx_Modules_b3.2.3.R" file and see line 1~20. Those are the list of R packages you may have to install unless you did. 

•	If you haven't installed any R pacakges below, it will take about 1 hour. 


## D. How to build iGenSig-Rx models based on CALGB dataset and predict the therapeutic response in ACOSOG and NOAH datasets.

* Find “iGenSigRx_RunCalculation_b3.2.3.R” file in your folder. This file contains the script to perform iGenSig-Rx modeling. 

* the iGenSig-Rx module has been tested on the latest R version 4.1.0 with Windows 10 computer Intel(R) Core(TM) i7-6700 CPU @ 3.40GHz, 32 GB RAM. 

* The iGenSig-Rx runs for one permutation and whole tranining set, instead of all 10 permutations. 

* Total running time on Windows 10 was ~ 2 hours. It can vary according to your computer's spec.

######################################################################################

>setwd("yourWorkingDirectory")

>source("iGenSigRx_Modules_b3.2.3.abbrv.R")

######################################################################################
##### Step1. Load treatment outcome data files of CALGB, ACOSOG, and NOAH.
######################################################################################
#phenoData must contain a column called "label" indicating whether the subject is "sen" or "res" or "mid"

> CALGB.phenotypefile<-"./PhenotypeData/CALGB_PhenotypesData_277s.tsv"
> CALGB.phenoData<-read.delim(CALGB.phenotypefile,stringsAsFactors=F,check.names=F,header=T,sep="\t")  
> CALGB.phenoData$SUBJECT_ID=gsub("(\\_)0+", "\\1", CALGB.phenoData$SUBJECT_ID, perl=TRUE)    
> CALGB.phenoData$label<-ifelse(CALGB.phenoData$pCR_Breast==1,"sen", ifelse(CALGB.phenoData$pCR_Breast==0,"res",NA))

> ACOSOG.phenotypefile<-"./PhenotypeData/ACOSOG_PhenotypeData_47s.tsv"   #ACOSOG_PhenotypeData.tsv"
> ACOSOG.phenoData<-read.delim(ACOSOG.phenotypefile, stringsAsFactors=F,check.names=F,header=T,sep="\t")
> ACOSOG.phenoData$label<-ifelse(ACOSOG.phenoData$pCR_Breast==1,"sen",ifelse(ACOSOG.phenoData$pCR_Breast==0,"res",NA)) #revise this based on different dataset

> NOAH.phenotypefile<-"./PhenotypeData/NOAH_PhenotypeData_63s.tsv" 
> NOAH.phenoData<-read.delim(NOAH.phenotypefile, stringsAsFactors=F,check.names=F,header=T,sep="\t")
> NOAH.phenoData$label<-ifelse(NOAH.phenoData$pCR_Breast==1,"sen",ifelse(NOAH.phenoData$pCR_Breast==0,"res",NA)) #revise this based on different dataset

#######################################################################################
##### Step2. Load CALGB genotype and preCal datasets for iGenSig analysis   
#######################################################################################
> TCGA.genotype.list <- readRDS("./GenotypePrecalmatrixData/TCGA.genotype.list.rds")
> CALGB.genotype.list <- readRDS("./GenotypePrecalmatrixData/TCGACALGBgenotype.list.rds")
> ACOSOG.genotype.list <- readRDS("./GenotypePrecalmatrixData/ACOSOG.genotype.list.rds")
> NOAH.genotype.list <- readRDS("./GenotypePrecalmatrixData/NOAH.genotype.list.rds")
> CALGB.preCalmatrix <- readRDS("./GenotypePrecalmatrixData/CALGB.preCalmatrix.rds")
> NOAH.preCalmatrix <- readRDS("./GenotypePrecalmatrixData/NOAH.preCalmatrix.rds")

#######################################################################################
#####  Step3. Load testset sample annotation data. 
#####		And Process GALGB.genotype and CALGB.phenoData to include only available subjects according to the annotation data
#######################################################################################
> FoldAssignFile="./TestsetAnnotationData/TestsetAssignments_CALGB_Arm1_2.tsv"
> fold.assign=read.delim(FoldAssignFile,stringsAsFactors=F,row.names=1, check.names=F, header=T, sep="\t")
> fold.assign["0",]=rep(0,ncol(fold.assign))# generate a set without testing set
> CALGB.genotype.list=lapply(CALGB.genotype.list,filter,filter=colnames(fold.assign),inlist=TRUE)  # 223,904
> CALGB.genotype.list=lapply(CALGB.genotype.list,function(x) x[x %in% colnames(fold.assign)]) 
> CALGB.phenoData=CALGB.phenoData[CALGB.phenoData$SUBJECT_ID %in% colnames(fold.assign),]  # 265  15
> CALGB.phenoData=CALGB.phenoData[!is.na(CALGB.phenoData$pCR_Breast),]

#######################################################################################
##### Step4. Parameter settings 
#######################################################################################
#please specify parameters: 

#1) the most important parameter is the sen.weightcut and res.weightcut. This parameter represents the signal to noise threshold.

#2) the power can be set between 1 or 2, and ECNpenalty can be set as 0.5 or 1. We used an universal parameter root=1, power=1, and ECNpenalty=1 in the manuscript.

> sen.weightcut<- 0.13; res.weightcut<- 0.13; root<- 1; ECNpenalty<- 1; power=1

#######################################################################################
##### Step5. Set up method and output directory for iGenSig-Rx
#######################################################################################
> parameters=list(genSig.outprefix="pearson",feature.select=NULL, confound=NULL,fold.assign=fold.assign,
>                 FoldAssignFile=FoldAssignFile,method="pearson",
>                 redun.method="ochiai",redun.cut=0.1,minsize=10,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,sen.pcut=1,res.pcut=1,
>                 power=power,root=root,rm.equivocal=TRUE,highlevel.minsize=15,ECNpenalty=ECNpenalty,by.cluster=FALSE,
>                 dgensig.method="math") 
> list2env(parameters, .GlobalEnv)
> OutTopDir<-"/zfs2/xiaosongwang/sal170/17_9_iGenSig_b3.0.3_NOAH_IJB_NeoALTTO/70_GithubSubmissionTest_iGenSigOnco/Result_iGenSigRx"
> dir.create(OutTopDir) 
> SubFolder<-paste0(OutTopDir, "/", sen.weightcut,"_RWC", res.weightcut, "_rt",root,"_ENCpn",ECNpenalty, "_Pw", power)
> dir.create(SubFolder)
> folder.prefix=paste0("iGenSig_CALGB")  
> CALGB.gensigdir=paste(SubFolder, "/",  folder.prefix,sep="")
> dir.create(CALGB.gensigdir)  # This is output sub-directory
> save(parameters,file=paste0(subject.gensigdir,"/iGenSig.parameters.rda"))

#######################################################################################
##### Step6. Calculate iGenSig-Rx score based on similarity method: CALGB clinical trial
#######################################################################################
#Perform weighted K-S tests for each permuted training/testing set or all CALGB subjects as training set
#If you want to calculate iGenSig-Rx scores for 10 permutations, please run the for look like "for (i in 1:nrow(fold.assign)) {

> for (i in c(1:11)) {			
>   outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
>   trainset=colnames(fold.assign)[which(fold.assign[i,]==0)]
>   trainset=trainset[trainset %in% CALGB.phenoData$SUBJECT_ID]
>   weightfile=paste0(subject.gensigdir,"/",outprefix,"_Weight.xls")
>   if (file.exists(weightfile)) {
>     print(paste("calcluating Gensig using",weightfile))
>     bachCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=trainset,outprefix=outprefix,
>                               subject.genotype.list=CALGB.genotype.list,preCalmatrix=CALGB.preCalmatrix,
>                               tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=CALGB.gensigdir,method=method,
>                               redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
>                               power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
>   }else{
>     subject.sensitive=CALGB.phenoData$SUBJECT_ID[(CALGB.phenoData$SUBJECT_ID %in% trainset) & CALGB.phenoData$pCR_Breast==1]
>     bachCal.GenSig.similarity(subject.sensitive=subject.sensitive,feature.select=feature.select,trainset=trainset,outprefix=outprefix,
>                               subject.genotype.list=CALGB.genotype.list,preCalmatrix=CALGB.preCalmatrix,
>                               tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=CALGB.gensigdir,method=method,
>                               redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
>                               power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,
>                               highlevel.minsize=highlevel.minsize,confound.factor=confound,by.cluster=by.cluster) #p.cut should be assigned NULL if similarity index is used instead of correlation statistics
>   }
> }

#######################################################################################
##### Step7. CALGB: dGenSig and benchmark test	 		
#######################################################################################
# We blocked this Step7 code lines which calculate iGenSig-RX performanc AUROC on CALGB because CALGB phenotypeData is not provided. CALGB phenotypeData is crendential.
#> batchCal.dGenSig.trial2trial.train (gensigdir = subject.gensigdir, test.matrix= fold.assign[,colnames(fold.assign) %in% CALGB.phenoData$SUBJECT_ID[CALGB.phenoData$Tx_arm!=3]],
#>   		phenoData = CALGB.phenoData[CALGB.phenoData$Tx_arm!=3,], method=dgensig.method)

#> batchbenchmark.dGenSig.trial(gensig.dir=subject.gensigdir, phenoData=CALGB.phenoData[CALGB.phenoData$Tx_arm!=3,], by.cluster=by.cluster)

#######################################################################################
##### Step8. preparation for ACOSOG analysis
#######################################################################################
> train.gensigdir=subject.gensigdir
> ACOSOG.gensigdir=paste0(train.gensigdir,"_ACOSOG")
> dir.create(ACOSOG.gensigdir)
> dgensig.model=mget(load(paste(train.gensigdir,"/iGenSig.parameters.rda",sep=""), envir=(tmp<- new.env())), envir=tmp)
> list2env(dgensig.model$parameters, .GlobalEnv)

#######################################################################################
##### Step9. Calculate iGenSig-Rx scores for Validation dataset: ACOSOG clinical trial data
#######################################################################################

#If you want to calculate iGenSig-Rx scores for 10 permutations, please run the for look like "for (i in 1:nrow(fold.assign)) {
> for (i in c(1:11)) {
>   outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
>   weightfile=paste0(train.gensigdir,"/",outprefix,"_Weight.xls")
>   print(paste("calcluating Gensig using",weightfile))
>   bachCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=NULL,outprefix=outprefix,
>                             subject.genotype.list=ACOSOG.genotype.list,preCalmatrix=ACOSOG.preCalmatrix,
>                             tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=ACOSOG.gensigdir,method=method,
>                             redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
>                             power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster,
>                             validationset=TRUE) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
> }
#######################################################################################
##### Step10. Calculate dGenSig scores and benchmark for Validation dataset: ACOSOG clinical trial data
#######################################################################################
> batchCal.dGenSig.trial2trial.validation (train.gensigdir=train.gensigdir, validation.gensigdir=ACOSOG.gensigdir,method=dgensig.method)
> batchbenchmark.dGenSig.trial(gensig.dir=ACOSOG.gensigdir, phenoData=ACOSOG.phenoData, by.cluster=by.cluster, event.col=15, time.col=16, validationset=TRUE)

#######################################################################################
##### Step11. preparation for NOAH analysis
#######################################################################################
> train.gensigdir=subject.gensigdir
> NOAH.gensigdir=paste0(train.gensigdir,"_NOAH")
> dir.create(NOAH.gensigdir)
> dgensig.model=mget(load(paste(train.gensigdir,"/iGenSig.parameters.rda",sep=""), envir=(tmp<- new.env())), envir=tmp)
> list2env(dgensig.model$parameters, .GlobalEnv)

#######################################################################################
##### Step12. Calculated GenSig scores for Validation data: NOAH clinical trial data
#######################################################################################
#If you want to calculate iGenSig-Rx scores for 10 permutations, please run the for look like "for (i in 1:nrow(fold.assign)) {

> for (i in c(1:11)) {
>   outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
>   weightfile=paste0(train.gensigdir,"/",outprefix,"_Weight.xls")
>   print(paste("calcluating Gensig using",weightfile))
>   bachCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=NULL,outprefix=outprefix,
>                             subject.genotype.list=NOAH.genotype.list,preCalmatrix=NOAH.preCalmatrix,
>                             tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=NOAH.gensigdir,method=method,
>                             redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
>                             power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster,
>                             validationset=TRUE) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
> }

#######################################################################################
##### Step12. option 1-2: calculate dGenSig and benchmark for Validation data: NOAH clinical trial data
#######################################################################################
> batchCal.dGenSig.trial2trial.validation (train.gensigdir=train.gensigdir, validation.gensigdir=NOAH.gensigdir, method=dgensig.method)
> batchbenchmark.dGenSig.trial(gensig.dir=NOAH.gensigdir, phenoData=NOAH.phenoData, by.cluster=by.cluster, validationset=TRUE)



## E. Expected outcome
•	Please find "./Result_iGenSigRx/0.13_RWC0.13_rt1_ENCpn1_Pw1/iGenSig_CALGB_ACOSOG/2022-07-02.dGenSig_XX_benchmark.result.xls" file. The file has the summary of the prediction performance AUROC.

