## DysRegSig
DysRegSig is capable of exploring gene dysregulations from high-dimensional data, ranking dysregulations and relevant TFs, and building explanatory signature based on gene dysregulations.

## Installation
DysRegSig dependents packages Boruta, limma, glmnet, ggpubr, ggplot2, igraph, expm, flare, RGBM, survival, and survcomp.
At first, a Bioconductor package needed to be pre-installed:

if(!require(BiocManager)) install.packages("BiocManager")

library(BiocManager)

install('survcomp')


Then use devtools to install DysRegSig package directly from github.

if(!require(devtools)) install.packages("devtools")

devtools::install_github('SCBIT-YYLab/DysRegSig')

## Input data


## User Guides


