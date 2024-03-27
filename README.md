---
title: "iGenSig-Rx"
author: "SanghoonLee"
date: "2024-03-26"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# How to use iGenSig-Rx R modules to build multi-omics models for predicting therapeutic responses based on GDSC dataset

## iGenSig-Rx: an integral genomic signature based white-box method for modeling clinical therapeutic responses using genome-wide sequencing data

• The current version of iGenSig-Rx is beta 3.2.3 (March 11th, 2024).

• The iGenSig-Rx was built on R version 4.1.0 and tested on Linux and Windows environment.

## A. Introduction

• Here we developed an integral genomic signature (iGenSig-Rx) analysis as a new class of transparent, interpretable, and resilient modeling methods for big data-based precision medicine with improved cross-dataset applicability and tolerance to sequencing bias.

• We define genomic features that significantly correlate with a clinical phenotype (such as therapeutic response) as genomic correlates, and an integral genomic signature as the integral set of redundant genomic correlates for a given clinical phenotype such as therapeutic response.

• We postulate that the redundant genomic features, which are usually eliminated through dimensionality reduction or feature removal during multi-omics modeling, may help overcome the sequencing bias.

• Using genomic dataset for chemical perturbations, we developed the iGenSig-Rx models predicting clinical response to HER2-targeted therapy and tested the cross-dataset performance of selected models based on independent multi-omics datasets for breast cancer HER2-positive patients assessing their therapeutic responses.

## B. How to build iGenSig-Rx models based on CALGB dataset and predict the therapeutic response in ACOSOG and NOAH datasets.

-   Find "iGenSigRx_RunCalculation_b3.2.3.R" file in your folder. This file contains the script to perform iGenSig-Rx modeling.

-   the iGenSig-Rx module has been tested on the latest R version 4.3.2 with MacOS Sonoma 14.4

-   The iGenSig-Rx runs for one permutation and whole tranining set, instead of all 10 permutations.

-   Total running time on Windows 10 was \~ 2 hours. It can vary according to your computer's spec.

## C. Installing R packages

```{Installation}
install.packages("devtools")
library(devtools)

devtools::install_github("wangxlab/iGenSig-Rx")
library(iGenSig-Rx)
```

## D. Running iGenSig-Rx

## Step1. Setting parameters and output directory

1)  the most important parameter is the sen.weightcut and res.weightcut. This parameter represents the signal to noise threshold. We have tested sen.weightcut and res.weightcut in the range 0.05\~2.5. We found that 0.13 was the best to obtain signal against noise.
2)  the power can be set between 1 or 2, and ECNpenalty can be set as 0.5 or 1. We used an universal parameter root=1, power=1, and ECNpenalty=1 in the manuscript.

```{Parameter setting}
sen.weightcut<- 0.13; res.weightcut<- 0.13; root<- 1; ECNpenalty<- 1; power=1

#Define the parameters dictionary
parameters=list(genSig.outprefix="pearson",feature.select=NULL, confound=NULL,fold.assign=fold.assign,
                 FoldAssignFile=FoldAssignFile,method="pearson",
                 redun.method="ochiai",redun.cut=0.1,minsize=10,sen.weightcut=sen.weightcut,
                 res.weightcut=res.weightcut,sen.pcut=1,res.pcut=1,power=power,root=root,
                 rm.equivocal=TRUE,highlevel.minsize=15,ECNpenalty=ECNpenalty,by.cluster=FALSE,
                 dgensig.method="math") 
list2env(parameters, .GlobalEnv)
OutTopDir<-"YourOutputDirectory"     ## <<<=========== You should set this up.
dir.create(OutTopDir) 
SubFolder<-paste0(OutTopDir,"/",sen.weightcut,"_RWC",res.weightcut,"_rt",root,"_ENCpn",ECNpenalty,"_Pw", power)
dir.create(SubFolder)
folder.prefix=paste0("iGenSig_CALGB")  
CALGB.gensigdir=paste(SubFolder, "/",  folder.prefix,sep="")
dir.create(CALGB.gensigdir)  # This is output sub-directory

#Save the parameters to a file
save(parameters,file=paste0(CALGB.gensigdir,"/iGenSig.parameters.rda"))
```

## Step2. Calculate iGenSig-Rx score based on similarity method: CALGB clinical trial

Perform weighted K-S tests for each permuted training/testing set or all CALGB subjects as training set

If you want to calculate iGenSig-Rx scores for 10 permutations, please run the for look like "for (i in 1:nrow(fold.assign)) {"

```{iGenSig_Rx score on CALGB}
for (i in c(1)) {
      outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
      trainset=colnames(fold.assign)[which(fold.assign[i,]==0)]
      trainset=trainset[trainset %in% CALGB.phenoData$SUBJECT_ID]
      weightfile=paste0(CALGB.gensigdir,"/",outprefix,"_Weight.xls")
      if (file.exists(weightfile)) {
          print(paste("calcluating Gensig using",weightfile))
          batchCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,
            trainset=trainset,outprefix=outprefix, subject.genotype.list=CALGB.genotype.list,
            preCalmatrix=CALGB.preCalmatrix,tcga.genotype.list=TCGA.genotype.list,
            subject.gensigdir=CALGB.gensigdir,method=method,redun.method=redun.method,redun.cut=redun.cut,
            minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,power=power,root=root,
            ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,highlevel.minsize=highlevel.minsize,
            by.cluster=by.cluster) 
            #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
      }else{
        subject.sensitive=CALGB.phenoData$SUBJECT_ID[(CALGB.phenoData$SUBJECT_ID %in% trainset) &
                                                       CALGB.phenoData$pCR_Breast==1]
        batchCal.GenSig.similarity(subject.sensitive=subject.sensitive,feature.select=feature.select,
            trainset=trainset,outprefix=outprefix,subject.genotype.list=CALGB.genotype.list,
            preCalmatrix=CALGB.preCalmatrix,tcga.genotype.list=TCGA.genotype.list,
            subject.gensigdir=CALGB.gensigdir,method=method,redun.method=redun.method,
            redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
            power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,
            highlevel.minsize=highlevel.minsize,confound.factor=confound,by.cluster=by.cluster) 
            #p.cut should be assigned NULL if similarity index is used instead of correlation statistics
  }
}
```

## Step3. CALGB: dGenSig and benchmark test

```{CALGB iGenSig_Rx score}
batchCal.dGenSig.trial2trial.train(gensigdir=CALGB.gensigdir,test.matrix=fold.assign[,colnames(fold.assign) %in%
       CALGB.phenoData$SUBJECT_ID],phenoData = CALGB.phenoData, method=dgensig.method)

batchbenchmark.dGenSig.trial(gensig.dir=subject.gensigdir, phenoData=CALGB.phenoData, by.cluster=by.cluster)
```

## Step4. Preparation for ACOSOG analysis

```{Prepare ACOSG analysis}
train.gensigdir=CALGB.gensigdir
ACOSOG.gensigdir=paste0(train.gensigdir,"_ACOSOG")
dir.create(ACOSOG.gensigdir)
dgensig.model=mget(load(paste(train.gensigdir,"/iGenSig.parameters.rda",sep=""), 
                        envir=(tmp<- new.env())), envir=tmp)
list2env(dgensig.model$parameters, .GlobalEnv)
```

## Step5. Option 1: calculate iGenSig-Rx scores for Validation dataset: ACOSOG clinical trial data \# this takes about 2 min.

```{calculate iGenSig_Rx score on ACOSG}
for (i in c(1)) {
    outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
    weightfile=paste0(train.gensigdir,"/",outprefix,"_Weight.xls")
    print(paste("calcluating Gensig using",weightfile))
    batchCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,
              trainset=NULL,outprefix=outprefix,subject.genotype.list=ACOSOG.genotype.list,
              preCalmatrix=ACOSOG.preCalmatrix,tcga.genotype.list=TCGA.genotype.list,
              subject.gensigdir=ACOSOG.gensigdir,method=method,redun.method=redun.method,
              redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
              power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,
              highlevel.minsize=highlevel.minsize,by.cluster=by.cluster,validationset=TRUE) 
              # p.cut should be assigned NULL if similarity index is used instead of correlation statistic
}
```

## Step6. Calculate dGenSig scores and benchmark for Validation dataset: ACOSOG clinical trial data

```{ACOSG_dGenSig}
batchCal.dGenSig.trial2trial.validation (train.gensigdir=train.gensigdir, validation.gensigdir=ACOSOG.gensigdir,
                 method=dgensig.method)
batchbenchmark.dGenSig.trial(gensig.dir=ACOSOG.gensigdir, phenoData=ACOSOG.phenoData, by.cluster=by.cluster, 
                event.col=15, time.col=16, validationset=TRUE)
```

## Step7. preparation for NOAH analysis

```{Prepare NOAH analysis}
train.gensigdir=CALGB.gensigdir
NOAH.gensigdir=paste0(train.gensigdir,"_NOAH")
dir.create(NOAH.gensigdir)
dgensig.model=mget(load(paste(train.gensigdir,"/iGenSig.parameters.rda",sep=""), 
                        envir=(tmp<- new.env())), envir=tmp)
list2env(dgensig.model$parameters, .GlobalEnv)
```

## Step8. Calculated GenSig scores for Validation data: NOAH clinical trial data

If you want to calculate iGenSig-Rx scores for 10 permutations, please run the for look like "for (i in 1:nrow(fold.assign)) {"

```{calculate iGenSig_Rx score on NOAH}
for (i in c(1)) {
   outprefix=paste0("trainset_", genSig.outprefix,"_Fold",rownames(fold.assign)[i])
   weightfile=paste0(train.gensigdir,"/",outprefix,"_Weight.xls")
   print(paste("calcluating Gensig using",weightfile))
   batchCal.GenSig.similarity(subject.sensitive=NULL,feature.select=feature.select,weightfile=weightfile,
            trainset=NULL,outprefix=outprefix,subject.genotype.list=NOAH.genotype.list,
            preCalmatrix=NOAH.preCalmatrix,tcga.genotype.list=TCGA.genotype.list,
            subject.gensigdir=NOAH.gensigdir,method=method,redun.method=redun.method,
            redun.cut=redun.cut,minsize=minsize,sen.weightcut=sen.weightcut,res.weightcut=res.weightcut,
            power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,
            highlevel.minsize=highlevel.minsize,by.cluster=by.cluster,validationset=TRUE) 
            #p.cut should be assigned NULL if similarity index is used instead of correlation statistic
 }
```

## Step9. calculate dGenSig and benchmark for Validation data: NOAH clinical trial data

```{NOAH_dGenSig}
batchCal.dGenSig.trial2trial.validation (train.gensigdir=train.gensigdir, 
                             validation.gensigdir=NOAH.gensigdir, method=dgensig.method)
batchbenchmark.dGenSig.trial(gensig.dir=NOAH.gensigdir, phenoData=NOAH.phenoData, 
                             by.cluster=by.cluster, validationset=TRUE)
```

## E. Expected outcome

Please find "./Result_iGenSigRx/iGenSig_CALGB_ACOSOG/2024-03-26.dGenSig_45_benchmark.result.xls" file. The file has the summary of the prediction performance AUROC.
