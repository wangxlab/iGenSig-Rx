# iGenSig-oncologist running takes 10~12 hours depending on your computer's specificity. 
# iGenSige wad devolped on R4.1.0
setwd("YourWorkingDirectory")
source("iGenSigOncologist_Modules_b3.2.3.R")

################################################################################################################   
###### Step1. Load treatment outcome data files of CALGB, ACOSOG, and NOAH.
################################################################################################################   
CALGB.phenotypefile<-"./PhenotypeData/CALGB_PhenotypesData_277s.tsv"
CALGB.phenoData<-read.delim(CALGB.phenotypefile,stringsAsFactors=F,check.names=F,header=T,sep="\t")  
CALGB.phenoData$SUBJECT_ID=gsub("(\\_)0+", "\\1", CALGB.phenoData$SUBJECT_ID, perl=TRUE)    
CALGB.phenoData$label<-ifelse(CALGB.phenoData$pCR_Breast==1,"sen", ifelse(CALGB.phenoData$pCR_Breast==0,"res",NA))

ACOSOG.phenotypefile<-"./PhenotypeData/ACOSOG_PhenotypeData_47s.tsv"   #ACOSOG_PhenotypeData.tsv"
ACOSOG.phenoData<-read.delim(ACOSOG.phenotypefile, stringsAsFactors=F,check.names=F,header=T,sep="\t")
ACOSOG.phenoData$label<-ifelse(ACOSOG.phenoData$pCR_Breast==1,"sen",ifelse(ACOSOG.phenoData$pCR_Breast==0,"res",NA)) #revise this based on different dataset

NOAH.phenotypefile<-"./PhenotypeData/NOAH_PhenotypeData_63s.tsv" 
NOAH.phenoData<-read.delim(NOAH.phenotypefile, stringsAsFactors=F,check.names=F,header=T,sep="\t")
NOAH.phenoData$label<-ifelse(NOAH.phenoData$pCR_Breast==1,"sen",ifelse(NOAH.phenoData$pCR_Breast==0,"res",NA)) #revise this based on different dataset

################################################################################################################
####### Step2. Load CALGB genotype and preCal datasets for iGenSig analysis   
################################################################################################################
TCGA.genotype.list <- readRDS("./GenotypePrecalmatrixData/TCGA.genotype.list.rds")
CALGB.genotype.list <- readRDS("./GenotypePrecalmatrixData/TCGACALGBgenotype.list.rds")
ACOSOG.genotype.list <- readRDS("./GenotypePrecalmatrixData/ACOSOG.genotype.list.rds")
NOAH.genotype.list <- readRDS("./GenotypePrecalmatrixData/NOAH.genotype.list.rds")
CALGB.preCalmatrix <- readRDS("./GenotypePrecalmatrixData/CALGB.preCalmatrix.rds")
NOAH.preCalmatrix <- readRDS("./GenotypePrecalmatrixData/NOAH.preCalmatrix.rds")

################################################################################################################
######  Step3. Load testset sample annotation data. 
######		And Process GALGB.genotype and CALGB.phenoData to include only available subjects according to the annotation data
################################################################################################################
FoldAssignFile="./TestsetAnnotationData/TestsetAssignments_CALGB_Arm1_2.tsv"
fold.assign=read.delim(FoldAssignFile,stringsAsFactors=F,row.names=1, check.names=F, header=T, sep="\t")
fold.assign["0",]=rep(0,ncol(fold.assign))# generate a set without testing set
CALGB.genotype.list=lapply(CALGB.genotype.list,filter,filter=colnames(fold.assign),inlist=TRUE)  # 223,904
CALGB.genotype.list=lapply(CALGB.genotype.list,function(x) x[x %in% colnames(fold.assign)]) 
CALGB.phenoData=CALGB.phenoData[CALGB.phenoData$SUBJECT_ID %in% colnames(fold.assign),]  # 265  15
CALGB.phenoData=CALGB.phenoData[!is.na(CALGB.phenoData$pCR_Breast),]

################################################################################################################
####### Step4. Parameter settings 
################################################################################################################
#please specify parameters: 
#1) the most important parameter is the sen.weightcut and res.weightcut. This parameter represents the signal to noise threshold.
#2) the power can be set between 1 or 2, and ECNpenalty can be set as 0.5 or 1. We used an universal parameter root=1, power=1, and ECNpenalty=1 in the manuscript.

sen.weightcut<- 0.13; res.weightcut<- 0.13; root<- 1; ECNpenalty<- 1; power=1

################################################################################################################
###### Step5. Set up method and output directory for iGenSig-oncologist
################################################################################################################
parameters=list(genSig.outprefix="pearson",feature.select=NULL, confound=NULL,fold.assign=fold.assign,
                FoldAssignFile=FoldAssignFile,method="pearson",
                redun.method="ochiai",redun.cut=0.1,minsize=10,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,sen.pcut=1,res.pcut=1,
                power=power,root=root,rm.equivocal=TRUE,highlevel.minsize=15,ECNpenalty=ECNpenalty,by.cluster=FALSE,
                dgensig.method="math") 
list2env(parameters, .GlobalEnv)
OutTopDir<-"/zfs2/xiaosongwang/sal170/17_9_iGenSig_b3.0.3_NOAH_IJB_NeoALTTO/70_GithubSubmissionTest_iGenSigOnco/Result_iGenSigOncologist"
dir.create(OutTopDir) 
folder.prefix=paste0("iGenSig_CALGB")  
CALGB.gensigdir=paste(OutTopDir, "/",  folder.prefix,sep="")
dir.create(CALGB.gensigdir)  # This is output sub-directory
save(parameters,file=paste0(subject.gensigdir,"/iGenSig.parameters.rda"))
################################################################################################################
####### Step5. Option 1 CALGB: calculate GenSig based on similarity method  		########################################
################################################################################################################

for (i in 1:2) {
  outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
  trainset=colnames(fold.assign)[which(fold.assign[i,]==0)]
  trainset=trainset[trainset %in% CALGB.phenoData$SUBJECT_ID]
  weightfile=paste0(subject.gensigdir,"/",outprefix,"_Weight.xls")
  if (file.exists(weightfile)) {
    print(paste("calcluating Gensig using",weightfile))
    bachCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=trainset,outprefix=outprefix,
                              subject.genotype.list=CALGB.genotype.list,preCalmatrix=CALGB.preCalmatrix,
                              tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=CALGB.gensigdir,method=method,
                              redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
                              power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
  }else{
    subject.sensitive=CALGB.phenoData$SUBJECT_ID[(CALGB.phenoData$SUBJECT_ID %in% trainset) & CALGB.phenoData$pCR_Breast==1]
    bachCal.GenSig.similarity(subject.sensitive=subject.sensitive,feature.select=feature.select,trainset=trainset,outprefix=outprefix,
                              subject.genotype.list=CALGB.genotype.list,preCalmatrix=CALGB.preCalmatrix,
                              tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=CALGB.gensigdir,method=method,
                              redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
                              power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,
                              highlevel.minsize=highlevel.minsize,confound.factor=confound,by.cluster=by.cluster) #p.cut should be assigned NULL if similarity index is used instead of correlation statistics
  }
}

################################################################################################################
####### Step6. CALGB: dGenSig and benchmark test	 		########################################
################################################################################################################
batchCal.dGenSig.trial2trial.train (gensigdir = subject.gensigdir, test.matrix= fold.assign[,colnames(fold.assign) %in% CALGB.phenoData$SUBJECT_ID[CALGB.phenoData$Tx_arm!=3]],
  		phenoData = CALGB.phenoData[CALGB.phenoData$Tx_arm!=3,], method=dgensig.method)

batchbenchmark.dGenSig.trial(gensig.dir=subject.gensigdir, phenoData=CALGB.phenoData[CALGB.phenoData$Tx_arm!=3,], by.cluster=by.cluster)

######################################################################################################
######## Step 7. preparation for ACOSOG analysis
######################################################################################################
train.gensigdir=subject.gensigdir
ACOSOG.gensigdir=paste0(train.gensigdir,"_ACOSOG")
dir.create(ACOSOG.gensigdir)
dgensig.model=mget(load(paste(train.gensigdir,"/iGenSig.parameters.rda",sep=""), envir=(tmp<- new.env())), envir=tmp)
list2env(dgensig.model$parameters, .GlobalEnv)
######################################################################################################
######## Step 8. Option 1 calculate GenSig scores for Validation dataset: ACOSOG clinical trial data
######################################################################################################
for (i in 1:2) {
  outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
  weightfile=paste0(train.gensigdir,"/",outprefix,"_Weight.xls")
  print(paste("calcluating Gensig using",weightfile))
  bachCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=NULL,outprefix=outprefix,
                            subject.genotype.list=ACOSOG.genotype.list,preCalmatrix=ACOSOG.preCalmatrix,
                            tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=ACOSOG.gensigdir,method=method,
                            redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
                            power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster,
                            validationset=TRUE) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
}
######################################################################################################
######## Step 9. Calculate dGenSig scores and benchmark for Validation dataset: ACOSOG clinical trial data
######################################################################################################
batchCal.dGenSig.trial2trial.validation (train.gensigdir=train.gensigdir, validation.gensigdir=ACOSOG.gensigdir,method=dgensig.method)

batchbenchmark.dGenSig.trial(gensig.dir=ACOSOG.gensigdir, phenoData=ACOSOG.phenoData, by.cluster=by.cluster, event.col=15, time.col=16, validationset=TRUE)

######################################################################################################
######## Step 10. preparation for NOAH analysis
######################################################################################################
train.gensigdir=subject.gensigdir
NOAH.gensigdir=paste0(train.gensigdir,"_NOAH")
dir.create(NOAH.gensigdir)
dgensig.model=mget(load(paste(train.gensigdir,"/iGenSig.parameters.rda",sep=""), envir=(tmp<- new.env())), envir=tmp)
list2env(dgensig.model$parameters, .GlobalEnv)

######################################################################################################
######## Step 11. option 1 calculated GenSig scores for Validation data: NOAH clinical trial data
######################################################################################################
for (i in 1:2) {
  outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
  weightfile=paste0(train.gensigdir,"/",outprefix,"_Weight.xls")
  print(paste("calcluating Gensig using",weightfile))
  bachCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=NULL,outprefix=outprefix,
                            subject.genotype.list=NOAH.genotype.list,preCalmatrix=NOAH.preCalmatrix,
                            tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=NOAH.gensigdir,method=method,
                            redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
                            power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster,
                            validationset=TRUE) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
}

######################################################################################################
######## Step 12. option 1-2: calculate dGenSig and benchmark for Validation data: NOAH clinical trial data
######################################################################################################
batchCal.dGenSig.trial2trial.validation (train.gensigdir=train.gensigdir,
                                         validation.gensigdir=NOAH.gensigdir,
                                         method=dgensig.method)

batchbenchmark.dGenSig.trial(gensig.dir=NOAH.gensigdir, phenoData=NOAH.phenoData, by.cluster=by.cluster, validationset=TRUE)

