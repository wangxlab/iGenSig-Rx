## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

NumbVector<-seq(1:10)
usethis::use_data(NumbVector, compress="xz", overwrite=TRUE)


CALGB.genotype.list <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/CALGB.genotype.list.rds")
usethis::use_data(CALGB.genotype.list, compress="xz", overwrite=TRUE)

CALGB.genotype.file <- "./data-raw/CALGB_GenotypeMatrixExample.gmt"
usethis::use_data(CALGB.genotype.file, compress="xz", overwrite=TRUE)


CALGB.genotype.list_Example <- read_gmt(CALGB.genotype.file, min=5)
usethis::use_data(CALGB.genotype.list_Example, compress="xz", overwrite=TRUE)

CALGB.genotype.matrix_Example <-genolist2matrix(CALGB.genotype.list_Example)
usethis::use_data(CALGB.genotype.matrix_Example, compress="xz", overwrite=TRUE)

################################################################################################################
###### Step1. Load clinical data files for clinical trials.
################################################################################################################
CALGB.phenotypefile<-"/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/PhenotypeData/CALGB_PhenotypesData_Original.tsv"
CALGB.phenoData<-read.delim(CALGB.phenotypefile,stringsAsFactors=F,check.names=F,header=T,sep="\t") %>% dplyr::select(c("dbGaP_Subject_ID","SUBJECT_ID","Tx_arm","pCR_Breast")) %>% dplyr::filter(Tx_arm!=3);
dim(CALGB.phenoData) # 217 4
colnames(CALGB.phenoData)[3:4] <-c("TreatType","TherapyRes")
CALGB.phenoData$SUBJECT_ID=gsub("(\\_)0+", "\\1", CALGB.phenoData$SUBJECT_ID, perl=TRUE)
CALGB.phenoData$label<-ifelse(CALGB.phenoData$TherapyRes==1,"sen", ifelse(CALGB.phenoData$TherapyRes==0,"res",NA))
## Don't run this;  usethis::use_data(CALGB.phenoData, compress="xz", overwrite=TRUE) #. I should process CALGB.phenodata again below.

ACOSOG.phenotypefile<-"/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/PhenotypeData/ACOSOG_PhenotypeData_47s.tsv"   #ACOSOG_PhenotypeData.tsv"
ACOSOG.phenoData<-read.delim(ACOSOG.phenotypefile, stringsAsFactors=F,check.names=F,header=T,sep="\t")
ACOSOG.phenoData$label<-ifelse(ACOSOG.phenoData$pCR_Breast==1,"sen",ifelse(ACOSOG.phenoData$pCR_Breast==0,"res",NA)) #revise this based on different dataset
usethis::use_data(ACOSOG.phenoData, compress="xz", overwrite=TRUE)

NOAH.phenotypefile<-"/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/PhenotypeData/NOAH_PhenotypeData_63s.tsv"
NOAH.phenoData<-read.delim(NOAH.phenotypefile, stringsAsFactors=F,check.names=F,header=T,sep="\t")
NOAH.phenoData$label<-ifelse(NOAH.phenoData$pCR_Breast==1,"sen",ifelse(NOAH.phenoData$pCR_Breast==0,"res",NA)) #revise this based on different dataset
usethis::use_data(NOAH.phenoData, compress="xz", overwrite=TRUE)

################################################################################################################
####### Step2. Load CALGB genotype and preCal datasets for iGenSig analysis
################################################################################################################
TCGA.genotype.list <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/TCGA.genotype.list.rds")
usethis::use_data(TCGA.genotype.list, compress="xz", overwrite=TRUE)

## CALGB.genotype.list is processed one more time.
CALGB.genotype.list <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/CALGB.genotype.list.rds")

ACOSOG.genotype.list <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/ACOSOG.genotype.list.rds")
usethis::use_data(ACOSOG.genotype.list, compress="xz", overwrite=TRUE)

NOAH.genotype.list <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/NOAH.genotype.list.rds")
usethis::use_data(NOAH.genotype.list, compress="xz", overwrite=TRUE)

CALGB.preCalmatrix <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/CALGB.preCalmatrix.rds")
usethis::use_data(CALGB.preCalmatrix, compress="xz", overwrite=TRUE)

ACOSOG.preCalmatrix <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/ACOSOG.preCalmatrix.rds")
usethis::use_data(ACOSOG.preCalmatrix, compress="xz", overwrite=TRUE)

NOAH.preCalmatrix <- readRDS("/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/GenotypePrecalmatrixData/NOAH.preCalmatrix.rds")
usethis::use_data(NOAH.preCalmatrix, compress="xz", overwrite=TRUE)

################################################################################################################
######  Step3. Load testset sample annotation data.
######		And Process GALGB.genotype and CALGB.phenoData to include only available subjects according to the annotation data
################################################################################################################
FoldAssignFile="/Volumes/Expansion/CCBR_XWANGLAB10_original/17_9_NOAH_IJB_NeoALTTO_ALTTO_ClinicalTrialData/70_2_GithubSubmissionActual_iGenSigRx/iGenSig-Rx/TestsetAnnotationData/TestsetAssignments_CALGB_Arm1_2.tsv"
fold.assign=read.delim(FoldAssignFile,stringsAsFactors=F,row.names=1, check.names=F, header=T, sep="\t")
fold.assign["0",]=rep(0,ncol(fold.assign))# generate a set without testing set
usethis::use_data(fold.assign, compress="xz", overwrite=TRUE)

CALGB.genotype.list=lapply(CALGB.genotype.list,FilterFunction,DataToFilter=colnames(fold.assign),inlist=TRUE)  # 223,904
CALGB.genotype.list=lapply(CALGB.genotype.list,function(x) x[x %in% colnames(fold.assign)])
usethis::use_data(CALGB.genotype.list, compress="xz", overwrite=TRUE)

subject.genotype.list <- CALGB.genotype.list
usethis::use_data(subject.genotype.list, compress="xz", overwrite=TRUE)

feature.list<-subject.genotype.list
usethis::use_data(feature.list, compress="xz", overwrite=TRUE)


CALGB.phenoData=CALGB.phenoData[CALGB.phenoData$SUBJECT_ID %in% colnames(fold.assign),]  # 265  15
#CALGB.phenoData=CALGB.phenoData[!is.na(CALGB.phenoData$pCR_Breast),] # This data is senstive. We can't open to public.

usethis::use_data(CALGB.phenoData, compress="xz", overwrite=TRUE)
usethis::use_data(DataToFilter, compress="xz", overwrite=TRUE)

################################################################################################################
####### Step4. Parameter settings
################################################################################################################
#please specify parameters:
#1) the most important parameter is the sen.weightcut and res.weightcut. This parameter represents the signal to noise threshold.
#2) the power can be set between 1 or 2, and ECNpenalty can be set as 0.5 or 1. We used an universal parameter root=1, power=1, and ECNpenalty=1 in the manuscript.

sen.weightcut<- 0.13; res.weightcut<- 0.13; root<- 1; ECNpenalty<- 1; power=1
usethis::use_data(sen.weightcut, compress="xz", overwrite=TRUE)
usethis::use_data(res.weightcut, compress="xz", overwrite=TRUE)
usethis::use_data(root, compress="xz", overwrite=TRUE)
usethis::use_data(ECNpenalty, compress="xz", overwrite=TRUE)
usethis::use_data(power, compress="xz", overwrite=TRUE)

################################################################################################################
###### Step5. Set up method and output directory for iGenSig-Rx
################################################################################################################
# parameters=c(genSig.outprefix="pearson",feature.select=NULL, confound=NULL,fold.assign=fold.assign,
#                 method="pearson",redun.method="ochiai",redun.cut=0.1,minsize=10,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,sen.pcut=1,res.pcut=1,
#                 power=power,root=root,rm.equivocal=TRUE,highlevel.minsize=15,ECNpenalty=ECNpenalty,by.cluster=FALSE,
#                 dgensig.method="math")
# usethis::use_data(parameters, compress="xz", overwrite=TRUE)
# list2env(parameters, .GlobalEnv)

genSig.outprefix="pearson";feature.select=NULL; confound=NULL; method="pearson";redun.method="ochiai";redun.cut=0.1;minsize=10;
sen.pcut=1;res.pcut=1;rm.equivocal=TRUE;highlevel.minsize=15;by.cluster=FALSE;dgensig.method="math"

usethis::use_data(genSig.outprefix, compress="xz", overwrite=TRUE); usethis::use_data(feature.select, compress="xz", overwrite=TRUE); usethis::use_data(confound, compress="xz", overwrite=TRUE)
usethis::use_data(method, compress="xz", overwrite=TRUE); usethis::use_data(redun.method, compress="xz", overwrite=TRUE); usethis::use_data(redun.cut, compress="xz", overwrite=TRUE)
usethis::use_data(minsize, compress="xz", overwrite=TRUE); usethis::use_data(sen.pcut, compress="xz", overwrite=TRUE); usethis::use_data(res.pcut, compress="xz", overwrite=TRUE);
usethis::use_data(rm.equivocal, compress="xz", overwrite=TRUE); usethis::use_data(highlevel.minsize, compress="xz", overwrite=TRUE); usethis::use_data(by.cluster, compress="xz", overwrite=TRUE); usethis::use_data(dgensig.method, compress="xz", overwrite=TRUE)

OutTopDir<-"./Result_iGenSigRx"
dir.create(OutTopDir)
usethis::use_data(OutTopDir, compress="xz", overwrite=TRUE)

folder.prefix=paste0("iGenSig_CALGB")
usethis::use_data(folder.prefix, compress="xz", overwrite=TRUE)

CALGB.gensigdir=paste(OutTopDir, "/",  folder.prefix,sep=""); print(CALGB.gensigdir)
dir.create(CALGB.gensigdir)  # This is output sub-directory
usethis::use_data(CALGB.gensigdir, compress="xz", overwrite=TRUE)


subject.gensigdir=CALGB.gensigdir
usethis::use_data(subject.gensigdir, compress="xz", overwrite=TRUE)


save(parameters,file=paste0(CALGB.gensigdir,"/iGenSig.parameters.rda"))

################################################################################################################
####### Step6. Option 1 CALGB: calculate GenSig based on similarity method  		## This takes about 15 min per permutation
################################################################################################################
for (i in c(1)) {
  outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
  trainset=colnames(fold.assign)[which(fold.assign[i,]==0)]
  trainset=trainset[trainset %in% CALGB.phenoData$SUBJECT_ID]
  weightfile=paste0(CALGB.gensigdir,"/",outprefix,"_Weight.xls")
  if (file.exists(weightfile)) {
    print(paste("calcluating Gensig using",weightfile))
    batchCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,trainset=trainset,outprefix=outprefix,
                               subject.genotype.list=CALGB.genotype.list,preCalmatrix=CALGB.preCalmatrix,
                               tcga.genotype.list=TCGA.genotype.list,subject.gensigdir=CALGB.gensigdir,method=method,
                               redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
                               power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,by.cluster=by.cluster) #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
  }else{
    subject.sensitive=CALGB.phenoData$SUBJECT_ID[(CALGB.phenoData$SUBJECT_ID %in% trainset) & CALGB.phenoData$pCR_Breast==1]
    batchCal.GenSig.similarity(subject.sensitive=subject.sensitive,feature.select=feature.select,trainset=trainset,outprefix=outprefix,
                               subject.genotype.list=CALGB.genotype.list,preCalmatrix=CALGB.preCalmatrix,
                               tcga.genotype.list=TCGA.genotype.list, subject.gensigdir=CALGB.gensigdir,method=method,
                               redun.method=redun.method,redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
                               power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,
                               highlevel.minsize=highlevel.minsize,confound.factor=confound,by.cluster=by.cluster) #p.cut should be assigned NULL if similarity index is used instead of correlation statistics
  }
}

usetihs::use_data(subject.sensitive, compress="xz", overwrite=TRUE)

list1<-subject.sensitive
usethis::use_data(list1, compress="xz", overwrite=TRUE)

list2<-CALGB.genotype.list
usethis::use_data(list2, compress="xz", overwrite=TRUE)

compare.list<-CALGB.genotype.list
usethis::use_data(compare.list, compress="xz", overwrite=TRUE)


list.all <- unique(unlist(CALGB.genotype.list))
usethis::use_data(list.all, compress="xz", overwrite=TRUE)

A=length(list1); B=length(list2)
tmp.weight=J/sqrt(A*B)
P=length(list.all)
J1=J-1;P1=P-1
tmp.weight.1=J1

usethis::use_data(tmp.weight, compress="xz", overwrite=TRUE)
usethis::use_data(tmp.weight.1, compress="xz", overwrite=TRUE)




usethis::use_data(outprefix, compress="xz", overwrite=TRUE)
usethis::use_data(trainset, compress="xz", overwrite=TRUE)
usethis::use_data(subject.sensitive, compress="xz", overwrite=TRUE)



# inlist=TRUE
# usethis::use_data(inlist, compress="xz", overwrite=TRUE)



