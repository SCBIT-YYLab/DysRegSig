## DysRegSig
DysRegSig is a robust tools to perform gene dysregulation analysis. The main improvement of DysRegSig is that it can consider the combinational effect of multiple regulators to target in gene dysregulation analysis. DysRegSig also offers tools to implement following analysis, including modules of exploring the association of individual gene dysregulaiton to clinical factors, ranking dysregulated regulations and TFs, building mechanistic signature with gene dysregulations by using genetic algorithm.

## Installation
DysRegSig dependents packages Boruta, RGBM, glmnet, expm, flare, limma, ggpubr, ggplot2, reshape2, igraph, survival, survcomp, survminer, ROCR, pROC, e1071, and so on. At first, Bioconductor packages need to be pre-installed.

```{r, eval = FALSE}
if(!require(BiocManager)) install.packages("BiocManager")
library(BiocManager)
install(c('limma','survcomp'))
```

Then use devtools to install DysRegSig package from github.

```{r, eval = FALSE}
if(!require(devtools)) install.packages("devtools")
devtools::install_github('SCBIT-YYLab/DysRegSig')
```

```{r, eval = FALSE}
library(DysRegSig)
```

## Quick start
In **DysRegSig**, the main function for gene dysregulation analysis is `DysReg`, which could idnetify gene dysregulations from high-dimensional data while considering the combinational effect of multiple regulators to target. DysReg first build conditional GRNs with tree-based feature selection algorithm, where  regulatory intensity and its confidential interval of each linkis estimated with a de-biased LASSO method. Gene dysregulations were then identified by integrating three properties including differential regulation, differential expression of target, and the consistency between differential regulation and differential expression.

```{r, eval = FALSE}
data(ExpData)
ExpData[1:5,1:5]

data(tf2tar)
head(tf2tar)

data(ClinData)
head(ClinData)

group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
exp.1 <- ExpData[,colnames(ExpData) %in% group.1]

group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
exp.2 <- ExpData[,colnames(ExpData) %in% group.2]

set.seed(1234)
test.genes <- sample(1:nrow(ExpData),1000)

tmp.1 <- exp.1[test.genes,]
tmp.2 <- exp.2[test.genes,]

dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
                     de.genes = NULL, de.pval = 0.05, 
                     grn.method = 'Boruta', 
                     pValue = 0.01, ci = 0.90, verbose = T)
                     
dysreg.res <- dysreg.out$dysreg
head(dysreg.res)

```

The expression pattern of two genes in one gene dysregulaiton could be visulazied by function `plotDysRegExp` (**Figure 1**). In order to more clearly visualize the differences of gene regulation between conditions, `plotDysRegExp` adds the regression lines and confidence interval shadows calculated by single variable regression for each condtion. 

```{r, eval = FALSE}
plotDysregExp(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1],
              exp.1 = exp.1,exp.2 = exp.2, 
              exp1.lab = 'Response',exp2.lab = 'No-response',
              dysreg = dysreg.res, method ='dysreg', conf.int.level = 0.95)
```

## Usage

This package offers offer several methods for gene dysreualtion analysis, and subsequent analysis toos for gene dysregulations. The more details about uasge of **DysRegSig** can be found: 

```{r, eval = FALSE}
vignette("DysRegSig")
```
