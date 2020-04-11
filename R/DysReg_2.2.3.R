
################################################################################
#' @name condiGRN
#' @title Build conditional GRN
#' @description Build conditional gene regulatory network (GRN) with expression data and prior reference network by using feature selection algorithm.
#' @usage 
#' condiGRN(exp.data, tf2tar, 
#'          method = 'Boruta', pValue = 0.01, 
#'          threshold = NULL, verbose = TRUE, ...)
#' 
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param tf2tar The prior reference GRN containing TF-target interaction pairs.
#' @param method The method used to build conditional GRN, such as 'Boruta', or 'RGBM'.
#' @param pValue Confidence level used in Boruta. Default value should be used.
#' @param threshold The threshould for weight in RGBM.
#' @param verbose A logical value indicating whether display the computating progress.
#' @param ... Other parameters passed to \link{Boruta} or \link{RGBM}.
#'  
#' @import Boruta
#' @import RGBM
#' @import igraph
#' @import reshape2
#' 
#' @details While using method Boruta, the predifined pValue could be used as the threshold to filter out nonsignificant regulatory relationships.  While using method RGBM, users need to explore the threshould of weight based on the output of RGBM to filter out nonsignificant regulatory relationships before following analysis.
#' 
#' @return Conditional GRN.
#' 
#' @references 
#' Kursa M B, Rudnicki W R. Feature Selection with the Boruta Package. Journal of Statistical Software. 2010, 36(11): 13.
#' @references 
#' Mall R, Cerulo L, Garofano L, et al. RGBM: regularized gradient boosting machines for identification of the transcriptional regulators of discrete glioma subtypes. Nucleic Acids Res. 2018, 46 (7), e39.
#' 
#' @examples
#' \donttest{
#' # Build a conditional GRN based on a reference GRN.
#' data(ExpData)
#' ExpData[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' set.seed(1234)
#' tmp <- exp.1[sample(1:nrow(exp.1),100),]
#' 
#' ## using method Boruta
#' net.1 <- condiGRN(exp.data = tmp, tf2tar = tf2tar, method = 'Boruta', pValue = 0.01)
#' 
#' }
#' 
#' @export
condiGRN <- function(exp.data, tf2tar, method = 'Boruta', pValue = 0.01, 
                     threshold = NULL, verbose = TRUE, ...){
  
  if (ncol(exp.data) < 5) {
    stop("Expression data must have at least five samples.")
  }
  if (ncol(exp.data) >= 5 & ncol(exp.data) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp.data[is.na(exp.data)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  tf2tar <- tf2tar[tf2tar$Target %in% rownames(exp.data),]
  tf2tar <- tf2tar[tf2tar$TF %in% rownames(exp.data),]
  
  ## build conditinal GRN with Boruta
  if (method == 'Boruta'){
    targets <- unique(tf2tar$Target)

    condition.GRN <- vector() 
    
    for(i in 1:length(targets)){
      reg.tar.i <- tf2tar[tf2tar$Target == targets[i],]
      reg.i <- unique(reg.tar.i$TF)
      
      if(length(reg.i) > 1){
        exp.i <- as.data.frame(t(exp.data[c(targets[i],reg.i),]))
        colnames(exp.i)[1] <- 'tar'
        
        boruta.res <- Boruta(tar~.,data = exp.i, pValue = pValue)
        dec.res <- data.frame(finical.des = boruta.res$finalDecision)
        dec.res$finical.des <- as.character(dec.res$finical.des)
        
        if(is.element('Tentative', dec.res$finical.des)){
          boruta.res <- TentativeRoughFix(boruta.res)
          reg.i <- getSelectedAttributes(boruta.res, withTentative = F)
        } else {
          reg.i <- rownames(dec.res)[which(dec.res$finical.des == 'Confirmed')]
        }
        
        reg.tar.i <- reg.tar.i[reg.tar.i$TF %in% reg.i,]
      }
      
      condition.GRN <- rbind(condition.GRN,reg.tar.i)
      
      
      if (verbose) {
        pb <- txtProgressBar(min = 0, max = length(targets), style = 3)
        setTxtProgressBar(pb,i)
      }
    }
  }
  
  ## build conditinal GRN with RGBM
  if (method == 'RGBM'){
   
    exp.data <- t(exp.data)
    
    net.matrix <- graph.data.frame(tf2tar)
    net.matrix <- as_adjacency_matrix(net.matrix)
    
    K <- matrix(0, nrow(exp.data), ncol(exp.data))
    colnames(K) <- colnames(exp.data)
    rownames(K) <- rownames(exp.data)
    
    condition.GRN <- RGBM(E = exp.data,K = K,
                          g_M = net.matrix, 
                          tfs = unique(tf2tar$TF),
                          targets = unique(tf2tar$Target))
    
    condition.GRN <- melt(condition.GRN, na.rm = TRUE)
    colnames(condition.GRN) <- c("TF", "Target", "weight")
    
    if(!is.null(weight)){
      condition.GRN <- condition.GRN[condition.GRN$weight > threshold, ]
    }
  }
  
  condition.GRN$TF <- as.character(condition.GRN$TF)
  condition.GRN$Target <- as.character(condition.GRN$Target)
  
  rownames(condition.GRN) <- NULL
  
  return(condition.GRN)
}



################################################################################
#' @name quantiReg
#' @title Quantify regulatory intensities of regulations
#' @description Quantify regulatory intensities and its confidence intervals with de-biased LASSO.
#' @usage 
#' quantiReg(exp.data, net, ci)
#' 
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param net Conditioanl GRN output from \code{\link{condiGRN}}.
#' @param ci  Confidence invetal of coefficient.
#' 
#' @import 
#' 
#' @return A data frame containing the regulatory intensity and its confidence invertal for each regulation.
#' 
#' @references 
#' Javanmard A, Montanari A. Confidence intervals and hypothesis testing for high-dimensional regression. Journal of Machine Learning Research. 2014, 15(1): 2869-909.
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' set.seed(1234)
#' tmp <- exp.1[sample(1:nrow(exp.1),100),]
#' 
#' ## using method Boruta
#' net.1 <- condiGRN(exp.data = tmp, tf2tar = tf2tar, method = 'Boruta', pValue = 0.01)
#' 
#' ## Quantify regulatory intensity
#' quanti.net.1 <- quantiReg(exp.data = tmp, net = net.1, ci = 0.90)
#' 
#' }
#' 
#' @export
quantiReg <- function(exp.data, net, ci = 0.95){
  
  #source('lasso_inference.r')
  
  if (ncol(exp.data) < 5) {
    stop("Expression data must have at least five samples.")
  }
  if (ncol(exp.data) >= 5 & ncol(exp.data) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp.data[is.na(exp.data)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  
  ## quantify regulatory intensity with de-baised LASSO
  net <- net[net$Target %in% rownames(exp.data),]
  net <- net[net$TF %in% rownames(exp.data),]
  
  targets <- unique(net$Target)
  
  quanti.reg <- do.call(rbind,lapply(1:length(targets),function(i,exp.data,net,ci){
    
    reg.tar.i <- net[net$Target == targets[i],]
    reg.i <- unique(reg.tar.i$TF)
    
    exp.i <- as.data.frame(t(exp.data[c(targets[i],reg.i),]),stringsAsFactors = F)
    colnames(exp.i)[1] <- 'tar'
    
    exp.i <- scale(exp.i,center = TRUE,scale = TRUE)
    
    if(ncol(exp.i)>5){
      
      lambda=vector()
      
      t <- 1
      while (t <= 5) {
        # first run glmnet
        lasso.fit <- cv.glmnet(x = exp.i[,-1],y = exp.i[,1],standardize = FALSE)
        lambda.i <- lasso.fit$lambda.min
        lambda <- c(lambda,lambda.i)
        
        t <- t+1
      }
      
      lambda <- median(lambda)
      
      lasso.fit.i <- SSLasso(X = exp.i[,-1],y = exp.i[,1],lambda = lambda,
                             alpha = 1-ci,verbose = F)
      
      ci.i <- data.frame(TF = names(lasso.fit.i$coef), Target = targets[i],
                         unb.coef = lasso.fit.i$unb.coef,
                         low.lim = lasso.fit.i$low.lim, up.lim = lasso.fit.i$up.lim,
                         P.val = lasso.fit.i$pvals, stringsAsFactors = F)
      
    } else {
      ff <- as.formula(paste0('tar ~ ',paste(colnames(exp.i)[-1], collapse = ' + ')))
      glm.fit <- glm(ff, data = as.data.frame(exp.i))
      glm.fit.coef <- summary(glm.fit)$coef
      glm.fit.ci <- suppressMessages(confint(glm.fit, level = ci,method = 'confint.glm'))
      
      ci.i <- data.frame(TF = rownames(glm.fit.coef)[-1], Target = targets[i],
                         unb.coef = glm.fit.coef[-1,1],
                         low.lim = glm.fit.ci[-1,1], up.lim = glm.fit.ci[-1,2],
                         P.val = glm.fit.coef[-1,4], stringsAsFactors = F)
      
    }
    
    return(ci.i)
    
  },exp.data = exp.data,net = net,ci = ci))
  
  rownames(quanti.reg) <- NULL
  
  return(quanti.reg)
  
}



################################################################################
#' @name DysReg
#' @title Identify gene dysregulations
#' @description Identify gene dysregulations by integrating three properties including differential regulation, differential expression of target, and the consistency between differential regulation and differential expression. DysReg could consider the combinational effect of multiple regulators to target expresion in gene dysregulation analysis.
#' @usage 
#' DysReg(exp.1, exp.2, tf2tar, 
#'        de.genes = NULL, de.pval = NULL, de.qval = NULL, de.logFC = NULL, 
#'        grn.method = 'Boruta', pValue = 0.01, threshold = NULL, 
#'        ci = 0.95, verbose = TRUE, ...)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param tf2tar The  prior reference GRN containing TF-target interaction pairs.
#' @param de.genes A dataframe for differential expression genes. If no de.genes offered, DysReg uses the default method limma to implement differential expression analysis. If de.genes offered, The dataframe must include three colomns, "GeneSymbol", "high.condition", "de.logFC". "high.condition" means which condition represents high expression level. "de.logFC" is the output logFC from differential expression analysis. 
#' @param de.pval The cutoff of pval for filtering differential expression genes. If you don't use this parameter to filter differential expression genes, this parameter could be set as NULL. If you use this parameter to filter differential expression genes, this parameter could be set as a special number, such as 0.05.
#' @param de.qval The cutoff of qval used for filtering differential expression genes.  If you don't use this parameter to filter differential expression genes, this parameter could be set as NULL. If you use this parameter to filter differential expression genes, this parameter could be set as a special number, such as 0.05.
#' @param de.logFC The cutoff of absolute logFC used for filtering differential expression genes. If you don't use this parameter to filter differential expression genes, this parameter could be set as NULL. If you use this parameter to filter differential expression genes, this parameter could be set as a special number, such as 0.5. This parameter could be used by combining with de.pval or de.qval.
#' @param grn.method The method used to build conditional GRN, such as 'Boruta', or 'RGBM'.
#' @param pValue Confidence level used in Boruta. Default value should be used.
#' @param threshold The threshould for weight in RGBM.
#' @param ci  The confidence invetal of coefficient.
#' @param verbose A logical value indicating whether display the computating progress.
#' @param ... Other parameters passed to \code{\link{condiGRN}}.
#' 
#' @details DysReg first builds conditional GRNs with tree-based feature selection algorithm, where regulatory intensity and its confidence interval of each link is estimated by a de-biased LASSO method. Gene dysregulations were then identified by integrating three properties including differential regulation, differential expression of target, and the consistency between differential regulation and differential expression.
#' 
#' @import limma
#' @import Boruta
#' @import RGBM
#' @import igraph
#' @import reshape2
#' 
#' @return 
#' The results of gene dysregulation analysis:
#'  \item{de.genes}{The differential expression genes between conditions.}
#'  \item{dysreg}{The identified gene dysregulations}
#'  
#' @seealso 
#'  \code{\link{condiGRN}}; \code{\link{quantiReg}}
#'  
#' @references 
#' Li Q, Li J, Dai W, et al. Differential regulation analysis reveals dysfunctional regulatory mechanism involving transcription factors and microRNAs in gastric carcinogenesis. Artif Intell Med. 2017, 77, 12-22.
#'  
#' @examples
#' \donttest{
#' data(ExpData)
#' ExpData[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' data(ClinData)
#' head(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', 
#'                      pValue = 0.01, ci = 0.90, verbose = T)
#' 
#' dysreg <- dysreg.out$dysreg
#' head(dysreg)
#' 
#' }
#' 
#' @export
DysReg <- function(exp.1, exp.2, tf2tar, de.genes = NULL,
                   de.pval = NULL, de.qval = NULL, de.logFC = NULL, 
                   grn.method = 'Boruta', pValue = 0.01, threshold = NULL, 
                   ci = 0.95,verbose = TRUE, ...){
  
  if (min(ncol(exp.1), ncol(exp.2)) < 5) {
    stop("each expression matrix must have at least five samples.")
  }
  
  if (ncol(exp.1) >= 5 & ncol(exp.1) < 10 | ncol(exp.2) >= 5 & ncol(exp.2) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp.1[is.na(exp.1)]) > 0 | length(exp.2[is.na(exp.2)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  if (class(exp.1) != "data.frame" | class(exp.2) != "data.frame") {
    exp.1 <- as.data.frame(exp.1)
    exp.2 <- as.data.frame(exp.2)
  }
  
  if(!is.null(de.genes)){
    need.cols <- c("GeneSymbol", "high.condition", "de.logFC")
    judged.res <- is.element(need.cols, colnames(de.genes))
    if("FALSE" %in% judged.res) {
      stop("The input of clin.data must include at least three columns, that are 'GeneSymbol', 'high.condition', 'de.logFC'. Please check the de.genes whether includes these columns and ensure the colname name is consistent with 'GeneSymbol', 'high.condition', 'de.logFC'")
    }
      
    warning("Please check the high expression condition in your de.genes data.")
  }
    
  if(is.null(de.genes)){
    if(verbose){
      cat('\n',paste('Identifying differential expression genes'),'\n')
    }
    
    ## filter differential expression genes
    exp.1$GeneSymbol <- rownames(exp.1)
    exp.2$GeneSymbol <- rownames(exp.2)
    merged.exp <- merge(exp.1, exp.2, by = 'GeneSymbol')
    rownames(merged.exp) <- merged.exp$GeneSymbol
    merged.exp <- merged.exp[,-1]
    
    exp.1 <- exp.1[,setdiff(colnames(exp.1),'GeneSymbol')]
    exp.2 <- exp.2[,setdiff(colnames(exp.2),'GeneSymbol')]
    
    data.group <- c(rep('condition1',dim(exp.1)[2]),rep('condition2',dim(exp.2)[2]))
    design <- model.matrix(~0+factor(data.group))
    colnames(design) <- levels(factor(data.group))
    rownames(design) <- colnames(merged.exp)
    contrast.matrix <- makeContrasts(paste0(c('condition1','condition2'),collapse = "-"),
                                     levels = design)
    
    fit <- lmFit(merged.exp,design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    diff.out <- topTable(fit2, coef = 1, number = nrow(merged.exp))
    diff.out <- na.omit(diff.out)
    
    diff.filter <- diff.out
    colnames(diff.filter)[c(1,4,5)] <- c('logFC','de.pval','de.qval')
    diff.filter$de.logFC <- abs(diff.filter$logFC)
    
    cutoff <- data.frame(de.logFC = NA,
                         de.pval = NA,
                         de.qval = NA, stringsAsFactors = F)
    
    cutoff$de.logFC <- de.logFC
    cutoff$de.pval <- de.pval
    cutoff$de.qval <- de.qval
    
    cutoff.s <- colnames(cutoff)
    
    if (is.element('de.logFC',cutoff.s)){
      diff.filter <- diff.filter[diff.filter$de.logFC > cutoff$de.logFC,]
      cutoff.s <- setdiff(colnames(cutoff), 'de.logFC')
      
      cutoff.type <- paste('param', length(cutoff.s), sep = '_')
      
      de.filter <- switch (
        cutoff.type,
        param_0 = diff.filter,
        param_1 = subset(diff.filter, 
                         diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,2]),
        param_2 = subset(diff.filter, 
                         diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,2] &
                           diff.filter[,match(cutoff.s[2],colnames(diff.filter))] < cutoff[,3]))
    } else {
      cutoff.type <- paste('param', length(cutoff.s), sep = '_')
      
      de.filter <- switch (
        cutoff.type,
        param_0 = diff.filter,
        param_1 = subset(diff.filter, 
                         diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,1]),
        param_2 = subset(diff.filter, 
                         diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,1] &
                           diff.filter[,match(cutoff.s[2],colnames(diff.filter))] < cutoff[,2]))
    }
    
    
    de.genes <- diff.out[rownames(diff.out) %in% rownames(de.filter),]
    de.genes$GeneSymbol <- rownames(de.genes)
    de.genes <- within(de.genes,{
      high.condition <- NA
      high.condition[logFC < 0] <- 2
      high.condition[logFC > 0] <- 1
    })
    
    de.genes <- de.genes[,c("GeneSymbol","high.condition","logFC")]
    colnames(de.genes)[3] <- 'de.logFC'
    rownames(de.genes) <- NULL
	
    rm(merged.exp)
    rm(data.group)
    rm(design)
    rm(contrast.matrix)
    rm(fit)
    rm(fit2)
    rm(diff.out)
    rm(cutoff)
    rm(cutoff.type)
    rm(cutoff.s)
    rm(diff.filter)
    rm(de.filter)
  }
  
  if(nrow(de.genes) < 1){
    stop('The number of DEGs is 0, please redefined cutoff for DEGs.')
    
  } else{
    
    if(verbose){
      cat('\n',paste('The number of DEGs:',nrow(de.genes)),'\n')
    }
  }
  
  if(verbose){
    cat('\n',paste('Building conditional GRN for exp.1'),'\n')
  }
  
  tf2tar.1 <- tf2tar[tf2tar$TF %in% rownames(exp.1),]
  tf2tar.1 <- tf2tar.1[tf2tar.1$Target %in% rownames(exp.1),]
  tf2tar.1 <- tf2tar.1[tf2tar.1$Target %in% de.genes$GeneSymbol,]
  
  net.1 <- condiGRN(exp.data = exp.1,tf2tar = tf2tar.1, method = grn.method, 
                    pValue = pValue, threshold = threshold)
  
  
  if(verbose){
    cat('\n',paste('Building conditional GRN for exp.2'),'\n')
  }

  tf2tar.2 <- tf2tar[tf2tar$TF %in% rownames(exp.2),]
  tf2tar.2 <- tf2tar.2[tf2tar.2$Target %in% rownames(exp.2),]
  tf2tar.2 <- tf2tar.2[tf2tar.2$Target %in% de.genes$GeneSymbol,]
  
  net.2 <- condiGRN(exp.data = exp.2,tf2tar = tf2tar.2, method = grn.method, 
                    pValue = pValue, threshold = threshold)
  
  net <- unique(rbind(net.1,net.2))
  
 
  if(verbose){
    cat('\n',paste('Quantifying regulatory intensity for exp.1'), '\n')
  }
  
  quanti.reg.1 <- quantiReg(exp.data = exp.1,net = net,ci = ci)
  
  if(verbose){
    cat('\n',paste('Quantifying regulatory intensity for exp.2'), '\n')
  }
  
  quanti.reg.2 <- quantiReg(exp.data = exp.2,net = net,ci = ci)
  
  
  dysreg.net <- merge(quanti.reg.1, quanti.reg.2, by = c('TF','Target'),all = T)
  
  colnames(dysreg.net)[3:10] <- c("unb.coef.1", "low.lim.1", "up.lim.1", "P.val.1", 
                                  "unb.coef.2", "low.lim.2", "up.lim.2", "P.val.2")

  dysreg.net[,3:10][is.na(dysreg.net[,3:10])] <- 0
  
  data.1 <- dysreg.net[dysreg.net$unb.coef.1 * dysreg.net$unb.coef.2 != 0,]
  data.2 <- dysreg.net[dysreg.net$unb.coef.1 * dysreg.net$unb.coef.2 == 0,]
  
  data.1 <- data.1[data.1$low.lim.1 > data.1$up.lim.2 | 
                     data.1$up.lim.1 < data.1$low.lim.2,]
  data.2 <- data.2[data.2$low.lim.1 * data.2$up.lim.1 > 0 | 
                     data.2$low.lim.2 * data.2$up.lim.2 > 0,]
  
  dysreg <- rbind(data.1,data.2)
  rm(data.1)
  rm(data.2)
  
  dysreg <- merge(dysreg,de.genes,by.x = 'Target',by.y = 'GeneSymbol')
  dysreg <- dysreg[,c(2,1,3:12)]
  
  dysreg <- rbind(
    subset(dysreg,
           dysreg$high.condition == 1 & (dysreg$unb.coef.2 - dysreg$unb.coef.1) < 0),
    subset(dysreg, 
           dysreg$high.condition == 2 & (dysreg$unb.coef.1 - dysreg$unb.coef.1) > 0))
  
  dysreg.res <- list(de.genes = de.genes,
                     dysreg = dysreg)
  
  if(verbose){
    cat('\n',paste('DysReg finished'),'\n')
  }
  
  return(dysreg.res)
}



################################################################################
#' @name DiffCor
#' @title Differential correlation analysis
#' @description Differential correlation analysis with Fishers' Z test.
#' @usage 
#' DiffCor(exp.1, exp.2, cor.method = 'pearson', p.adj)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param cor.method The method used to calculate correlation coefficient.
#' @param p.adj  Correction method for p value adjust.
#' 
#' @return The identified gene pairs with significantly differential correlation.
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' ## implement differential correlation analysis
#' diffcor.res <- DiffCor(exp.1 = tmp.1, exp.2 = tmp.2, 
#'                        cor.method = 'pearson', p.adj = 'BH')
#'                        
#' ## set cutoff
#' diffcor.res <- subset(diffcor.res,p.val < 0.01)
#' head(diffcor.res)
#' 
#' }
#' 
#' @export
DiffCor <- function(exp.1, exp.2, 
                    cor.method = 'pearson', p.adj = 'BH'){
  
  if (min(ncol(exp.1), ncol(exp.2)) < 5) {
    stop("each expression matrix must have at least five samples.")
  }
  
  if (ncol(exp.1) >= 5 & ncol(exp.1) < 10 | ncol(exp.2) >= 5 & ncol(exp.2) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp.1[is.na(exp.1)]) > 0 | length(exp.2[is.na(exp.2)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  if (class(exp.1) != "matrix" | class(exp.2) != "matrix") {
    exp.1 <- as.matrix(exp.1)
    exp.2 <- as.matrix(exp.2)
  }
  
  merged.genes <- intersect(rownames(exp.1),rownames(exp.2))

  if(length(rownames(exp.1)) != length(merged.genes)| 
     length(rownames(exp.2)) != length(merged.genes)){
    warning("The gene list in exp.1 and exp.2 are not same.")
  }
  
  data.1 <- data.frame(genes = rownames(exp.1), exp.1)
  data.2 <- data.frame(genes = rownames(exp.2), exp.2)
  data <- merge(data.1, data.2, by = 'genes')
  rownames(data) <- data$genes
  
  rm(data.1)
  rm(data.2)
  
  exp.1 <- data[ , colnames(exp.1)]
  exp.2 <- data[ , colnames(exp.2)]
  
  cor.1 <- cor(t(exp.1), method = cor.method)
  cor.2 <- cor(t(exp.2), method = cor.method)
  
  z.1 <- (0.5*log((1 + cor.1)/(1 - cor.1)))
  z.2 <- (0.5*log((1 + cor.2)/(1 - cor.2)))
  
  z.value <- (z.2 - z.1)/((1/(ncol(exp.1)-3) + 1/(ncol(exp.2)-3))^0.5)
  
  gene.list <- as.character(rownames(data))
  
  diffcor <- do.call(rbind,lapply(1 : (length(gene.list) - 1), 
                               function(i, gene.list, cor.1, cor.2, z.value){
    
    cor1.i <- cor.1[rownames(cor.1) == gene.list[i], (i+1):length(gene.list)]
    cor1.i <- as.data.frame(cor1.i)
    cor1.i$gene.1 <- gene.list[i]
    cor1.i$gene.2 <- rownames(cor1.i)
    
    cor2.i <- cor.2[rownames(cor.2) == gene.list[i], (i+1):length(gene.list)]
    cor2.i <- as.data.frame(cor2.i)
    cor2.i$gene.1 <- gene.list[i]
    cor2.i$gene.2 <- rownames(cor2.i)
    
    diffcor.i <- merge(cor1.i,cor2.i,by = c('gene.1','gene.2'))
    
    z.i <- z.value[rownames(z.value)==gene.list[i], (i+1):length(gene.list)]
    z.i <- as.data.frame(z.i)
    z.i$gene.1 <- gene.list[i]
    z.i$gene.2 <- rownames(z.i)
    
    diffcor.i <- merge(diffcor.i,z.i,by = c('gene.1','gene.2'))
    
    diffcor.i <- diffcor.i[diffcor.i$gene.1 != diffcor.i$gene.2,]
    
    colnames(diffcor.i)[3:5] <- c('cor.1','cor.2','z.value')
    rownames(diffcor.i) <- NULL
    
    return(diffcor.i)
    
  }, gene.list = gene.list, cor.1 = cor.1, cor.2=cor.2, z.value=z.value))
  
  diffcor$p.val <- 2 * (1 - pnorm(abs(diffcor$z.value)))
  diffcor$p.adj <- p.adjust(diffcor$p.val,method = p.adj)
  
  diffcor <- diffcor[order(abs(diffcor$z.value),decreasing = T),]
  rownames(diffcor) <- NULL
  
  return(diffcor) 
}



################################################################################
#' @name DiffReg
#' @title Differential regulation analysis
#' @description Differential regulation analysis with Fishers' Z test.
#' @usage 
#' DiffReg(exp.1, exp.2, tf2tar, cor.method = 'pearson', p.adj)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param tf2tar The  prior reference GRN containing TF-target interaction pairs.
#' @param cor.method The method used to calculate correlation coefficient.
#' @param p.adj  Correction method for p value adjust.
#' 
#' @return The identified differential regulations.
#' 
#' @examples
#' \donttest{
#' 
#' data(ExpData) 
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' ## implement differential regulation analysis
#' diffreg.res <- DiffReg(tmp.1, tmp.2, tf2tar, 
#'                        cor.method = 'pearson', p.adj = 'BH')
#'                        
#' ## set cutoff
#' diffreg.res <- subset(diffreg.res,p.val < 0.01)
#' head(diffreg.res)
#' 
#' }
#' 
#' @export
DiffReg <- function(exp.1, exp.2, tf2tar, 
                    cor.method = 'pearson', p.adj = 'BH'){
  
  if (min(ncol(exp.1), ncol(exp.2)) < 5) {
    stop("each expression matrix must have at least five samples.")
  }
  
  if (ncol(exp.1) >= 5 & ncol(exp.1) < 10 | ncol(exp.2) >= 5 & ncol(exp.2) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp.1[is.na(exp.1)]) > 0 | length(exp.2[is.na(exp.2)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }

  if (class(exp.1) != "matrix" | class(exp.2) != "matrix") {
    exp.1 <- as.matrix(exp.1)
    exp.2 <- as.matrix(exp.2)
  }
  
  merged.genes <- intersect(rownames(exp.1),rownames(exp.2))
  
  if(length(rownames(exp.1)) != length(merged.genes)| 
     length(rownames(exp.2)) != length(merged.genes)){
    warning("exp.1 and exp.2 contain different gene list.")
    
    exp.1 <- exp.1[rownames(exp.1) %in% merged.genes,]
    exp.2 <- exp.2[rownames(exp.2) %in% merged.genes,]
  }
  
  exp.1 <- exp.1[order(rownames(exp.1), decreasing = F),]
  exp.2 <- exp.2[order(rownames(exp.2), decreasing = F),]
  

  tf2tar <- tf2tar[tf2tar$TF %in% merged.genes,]
  tf2tar <- tf2tar[tf2tar$Target %in% merged.genes,]
  
  tf.list <- unique(tf2tar$TF)
  
  cor.1 <- cor(t(exp.1),method = cor.method)
  cor.1 <- cor.1[rownames(cor.1) %in% tf.list,]
  
  cor.2 <- cor(t(exp.2),method = cor.method)
  cor.2 <- cor.2[rownames(cor.2) %in% tf.list,]
  
  z.1 <- (0.5*log((1 + cor.1)/(1 - cor.1)))
  z.2 <- (0.5*log((1 + cor.2)/(1 - cor.2)))
  
  z.value <- (z.2 - z.1)/((1/(ncol(exp.1)-3) + 1/(ncol(exp.2)-3))^0.5)
  
  diffreg <- do.call(rbind,lapply(1:length(tf.list), 
                                  function(i, tf2tar, tf.list, cor.1, cor.2, z.value){
    
    cor1.i <- cor.1[rownames(cor.1) == tf.list[i],]
    cor1.i <- as.data.frame(cor1.i)
    cor1.i$TF=tf.list[i]
    cor1.i$Target=rownames(cor1.i)
    
    cor2.i <- cor.2[rownames(cor.2) == tf.list[i],]
    cor2.i <- as.data.frame(cor2.i)
    cor2.i$TF <- tf.list[i]
    cor2.i$Target <- rownames(cor2.i)
    
    diffreg.i <- merge(cor1.i,cor2.i,by = c('TF','Target'))
    
    z.i <- z.value[rownames(z.value)==tf.list[i],]
    z.i <- as.data.frame(z.i)
    z.i$TF <- tf.list[i]
    z.i$Target <- rownames(z.i)
    
    diffreg.i <- merge(diffreg.i,z.i,by = c('TF','Target'))
    
    tf.i <- tf2tar[tf2tar$TF == tf.list[i],]
    diffreg.i <- diffreg.i[diffreg.i$Target %in% tf.i$Target,]
    colnames(diffreg.i)[3:5] <- c('cor.1','cor.2','z.value')
    rownames(diffreg.i) <- NULL
    
    return(diffreg.i)
    
  }, tf2tar = tf2tar,tf.list = tf.list,cor.1 = cor.1, cor.2 = cor.2, z.value = z.value))
  
  diffreg$p.val <- 2 * (1 - pnorm(abs(diffreg$z.value)))
  diffreg$p.adj <- p.adjust(diffreg$p.val,method = p.adj)
  
  diffreg <- diffreg[order(abs(diffreg$z.value),decreasing = T),]
  rownames(diffreg) <- NULL
  
  return(diffreg) 
}



################################################################################
#' @name DiffRegPlus
#' @title Implement differential regulation analysis based on \code{\link{DiffReg}} by combining differential expression of target, and the consistency between correlation change and target expression change
#' 
#' @description Identify regulations with differential correlation, differential expression of target, and the consistency between differential correlation and differential expression.
#' 
#' @usage 
#' DiffRegPlus(exp.1, exp.2, tf2tar, 
#'             de.genes = NULL, de.pval = NULL, de.qval = NULL, de.logFC = NULL, 
#'             cor.method = "pearson", p.adj = "BH", verbose = TRUE)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' 
#' @param tf2tar The  prior reference GRN containing TF-target interaction pairs.
#' @param de.genes A dataframe for differential expression genes. If no de.genes offered, DysReg uses the default method limma to implement differential expression analysis. If de.genes offered, The dataframe must include three colomns, "GeneSymbol", "high.condition", "de.logFC". "high.condition" means which condition represents high expression level. "de.logFC" is the output logFC from differential expression analysis. 
#' @param de.pval The cutoff of pval for filtering differential expression genes. If you don't use this parameter to filter differential expression genes, this parameter could be set as NULL. If you use this parameter to filter differential expression genes, this parameter could be set as a special number, such as 0.05.
#' @param de.qval The cutoff of qval used for filtering differential expression genes.  If you don't use this parameter to filter differential expression genes, this parameter could be set as NULL. If you use this parameter to filter differential expression genes, this parameter could be set as a special number, such as 0.05.
#' @param de.logFC The cutoff of absolute logFC used for filtering differential expression genes. If you don't use this parameter to filter differential expression genes, this parameter could be set as NULL. If you use this parameter to filter differential expression genes, this parameter could be set as a special number, such as 0.5. This parameter could be used by combining with de.pval or de.qval.
#' @param cor.method Which correlation coefficient (or covariance) is to be computed. One of "pearson" (default) or "spearman", can be abbreviated.
#' @param p.adj  Correction method for p value adjust.
#' @param verbose A logical value indicating whether display the computing progress.
#' 
#' @import limma
#' 
#' @return The identified regulations.
#' 
#' @seealso 
#' \code{\link{DiffReg}}
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' ## implement differential regulation analysis
#' diffreg.p.res <- DiffRegPlus(tmp.1,tmp.2, tf2tar, 
#'                              de.genes = NULL, de.pval = 0.05, 
#'                              cor.method = 'pearson', p.adj = 'BH')
#' 
#' ## set cutoff
#' diffreg.p.res <- subset(diffreg.p.res,p.val < 0.05)
#' head(diffreg.p.res)
#' 
#' }
#' 
#' @export
DiffRegPlus <- function(exp.1,exp.2, tf2tar,  
                        de.genes = NULL, 
                        de.pval=NULL, de.qval = NULL, de.logFC = NULL,
                        cor.method = 'pearson', p.adj = 'BH', verbose = TRUE){
  
  if (min(ncol(exp.1), ncol(exp.2)) < 5) {
    stop("each expression matrix must have at least five samples.")
  }
  
  if (ncol(exp.1) >= 5 & ncol(exp.1) < 10 | ncol(exp.2) >= 5 & ncol(exp.2) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp.1[is.na(exp.1)]) > 0 | length(exp.2[is.na(exp.2)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  if (class(exp.1) != "data.frame" | class(exp.2) != "data.frame") {
    exp.1 <- as.data.frame(exp.1)
    exp.2 <- as.data.frame(exp.2)
  }
  
  if(!is.null(de.genes)){
    need.cols <- c("GeneSymbol", "high.condition", "de.logFC")
    judged.res <- is.element(need.cols, colnames(de.genes))
    if("FALSE" %in% judged.res) {
      stop("The input of clin.data must include at least three columns, that are 'GeneSymbol', 'high.condition', 'de.logFC'. Please check the de.genes whether includes these columns and ensure the colname name is consistent with 'GeneSymbol', 'high.condition', 'de.logFC'")
    }
    
    warning("Please check the high expression condition in your de.genes data.")
  }
  
  if(is.null(de.genes)){
    if(verbose){
      cat('\n',paste('Identifying differential expression genes'),'\n')
    }
    
    ## filter differential expression genes
    exp.1$GeneSymbol <- rownames(exp.1)
    exp.2$GeneSymbol <- rownames(exp.2)
    merged.exp <- merge(exp.1, exp.2, by = 'GeneSymbol')
    rownames(merged.exp) <- merged.exp$GeneSymbol
    merged.exp <- merged.exp[,-1]
    
    exp.1 <- exp.1[,setdiff(colnames(exp.1),'GeneSymbol')]
    exp.2 <- exp.2[,setdiff(colnames(exp.2),'GeneSymbol')]
    
    data.group <- c(rep('condition1',dim(exp.1)[2]),rep('condition2',dim(exp.2)[2]))
    design <- model.matrix(~0+factor(data.group))
    colnames(design) <- levels(factor(data.group))
    rownames(design) <- colnames(merged.exp)
    contrast.matrix <- makeContrasts(paste0(c('condition1','condition2'),collapse = "-"),
                                     levels = design)
    
    fit <- lmFit(merged.exp,design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    diff.out <- topTable(fit2, coef = 1, number = nrow(merged.exp))
    diff.out <- na.omit(diff.out)
    
    diff.filter <- diff.out
    colnames(diff.filter)[c(1,4,5)] <- c('logFC','de.pval','de.qval')
    diff.filter$de.logFC <- abs(diff.filter$de.logFC)
    
    cutoff <- data.frame(de.logFC = NA,
                         de.pval = NA,
                         de.qval = NA, stringsAsFactors = F)
    
    cutoff$de.logFC <- de.logFC
    cutoff$de.pval <- de.pval
    cutoff$de.qval <- de.qval
    
    cutoff.s <- colnames(cutoff)
    
    if (is.element('de.logFC',cutoff.s)){
      diff.filter <- diff.filter[diff.filter$de.logFC > cutoff$de.logFC,]
      cutoff.s <- setdiff(colnames(cutoff), 'de.logFC')
      
      cutoff.type <- paste('param', length(cutoff.s), sep = '_')
      
      de.filter <- switch (
	    cutoff.type,
        param_0 = diff.filter,
        param_1 = subset(diff.filter, 
		                 diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,2]),
        param_2 = subset(diff.filter, 
		                 diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,2] & 
                           diff.filter[,match(cutoff.s[2],colnames(diff.filter))] < cutoff[,3]))
						   
    } else {
      cutoff.type <- paste('param', length(cutoff.s), sep = '_')
      
      de.filter <- switch (
	    cutoff.type,
        param_0 = diff.filter,
        param_1 = subset(diff.filter, 
		                 diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,1]),
        param_2 = subset(diff.filter, 
		                 diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,1] & 
                         diff.filter[,match(cutoff.s[2],colnames(diff.filter))] < cutoff[,2]))
    }
    
    de.genes <- diff.out[rownames(diff.out) %in% rownames(de.filter),]
    de.genes$GeneSymbol <- rownames(de.genes)
    de.genes <- within(de.genes,{
      high.condition <- NA
      high.condition[logFC < 0] <- 2
      high.condition[logFC > 0] <- 1
    })
	
    de.genes <- de.genes[,c("GeneSymbol","high.condition","logFC")]
    colnames(de.genes)[3] <- 'de.logFC'
    rownames(de.genes) <- NULL
	
    rm(merged.exp)
    rm(data.group)
    rm(design)
    rm(contrast.matrix)
    rm(fit)
    rm(fit2)
    rm(diff.out)
    rm(cutoff)
    rm(cutoff.type)
    rm(cutoff.s)
    rm(diff.filter)
    rm(de.filter)
  }
  
  if(nrow(de.genes) < 1){
    stop('The number of DEGs is 0, please redefined cutoff for DEGs.')
  } else{
    if(verbose){
      cat('\n',paste('The number of DEGs:',nrow(de.genes)),'\n')
	}
  }
  
  if(verbose){
    cat('\n',paste('Implementing differential regulation analysis'),'\n')
  }
  
  diffreg.res <- DiffReg(exp.1, exp.2, tf2tar, cor.method = cor.method, p.adj = p.adj)

  diffreg.plus <- merge(diffreg.res,de.genes,by.x = 'Target',by.y = 'GeneSymbol')
  
  diffreg.plus <- diffreg.plus[,c(2,1,3:9)]
  diffreg.plus <- rbind(
                    subset(diffreg.plus, 
					       diffreg.plus$high.condition == 1 & diffreg.plus$z.value < 0),
                    subset(diffreg.plus, 
					       diffreg.plus$high.condition == 2 & diffreg.plus$z.value > 0))

  diffreg.plus <- diffreg.plus[order(abs(diffreg.plus$z.value),decreasing = T),]
  rownames(diffreg.plus) <- NULL
  
  return(diffreg.plus)
}



###############################################################################
#' @name plotDysregExp
#' @title Visualize expression pattern of a gene dysreglaiton between conditions
#' @description Visualize expression pattern of a gene dysreglaiton between conditions.
#' @usage  
#' plotDysregExp(tf, tar, exp.1, exp.2, 
#'               exp1.label, exp2.label, 
#'               method, dysreg, conf.int.level, 
#'               cor.method = 'pearson', ...)
#' 
#' @param tf The TF of a gene dysregulation.
#' @param tar The target of a gene dysregulation.
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp1.label The label condition of exp.1.
#' @param exp2.label The label condition of exp.2.
#' @param method The method of dysregulation analysis, "DysReg", "DiffCor", "DiffReg", or "DiffRegPlus".
#' @param dysreg If method is "DysReg", this parameter offers the results of \code{\link{DysReg}}.
#' @param conf.int.level Level controlling confidence region. 
#' @param cor.method Which correlation coefficient (or covariance) is to be computed. One of "pearson" (default) or "spearman", can be abbreviated. the parameter is used while parameter method is "DiffCor", "DiffReg", or "DiffRegPlus".
#' @param ... Other parameters passed to \link{ggscatter}.
#' 
#' @import ggpubr
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', pValue = 0.01, ci = 0.90)
#' 
#' dysreg.res <- dysreg.out$dysreg
#' 
#' plotDysregExp(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1],
#'               exp.1 = exp.1, exp.2 = exp.2, 
#'               exp1.label = 'Response', exp2.label = 'No-response',
#'               dysreg = dysreg.res, method ='DysReg', 
#'               conf.int.level = 0.95)
#'               
#' }
#' 
#' @export
plotDysregExp=function (tf, tar, exp.1, exp.2, exp1.label, exp2.label, method, 
                        dysreg, conf.int.level, cor.method = 'pearson', ...) {
  
  data = rbind(data.frame(TF = as.numeric(exp.1[rownames(exp.1) == tf, ]), 
                          Target = as.numeric(exp.1[rownames(exp.1) == tar, ]), 
                          Group = exp1.label, stringsAsFactors = F), 
               data.frame(TF = as.numeric(exp.2[rownames(exp.2) == tf, ]), 
                          Target = as.numeric(exp.2[rownames(exp.2) == tar, ]), 
                          Group = exp2.label, stringsAsFactors = F))
  
  data$Group <- factor(data$Group,levels = c(exp1.label,exp2.label))
  
  if (method == "DysReg") {
    reg.i = dysreg[dysreg$TF == tf & dysreg$Target == tar, ]
    
    sp <- ggscatter(data, x = "TF", y = "Target", add = "reg.line", 
                    conf.int = TRUE, conf.int.level = conf.int.level, 
                    color = "Group", shape = "Group", 
                    palette = rainbow(5)[c(1, 4)]) + 
      xlab(tf) + ylab(tar)
  }
  else {
    sp <- ggscatter(data, x = "TF", y = "Target", add = "reg.line", 
                    conf.int = TRUE, color = "Group", shape = "Group", 
                    palette = rainbow(5)[c(1, 4)]) + 
    xlab(tf) + ylab(tar)
  }
  
  sp <- switch(method, 
               DiffCor = sp + 
                 stat_cor(aes(color = Group), 
                          method = cor.method, size = 6, 
                          label.x = min(data$TF) + (max(data$TF) - min(data$TF)) * 0.05, 
                          label.y = c(max(data$Target) + 
                                        (max(data$Target) - min(data$Target)) * 0.2, 
                                      max(data$Target) + 
                                        (max(data$Target) - min(data$Target)) * 0.1)), 
               
               DiffReg = sp + 
                 stat_cor(aes(color = Group), 
                          method = cor.method, size = 6, 
                          label.x = min(data$TF) + (max(data$TF) - min(data$TF)) * 0.05, 
                          label.y = c(max(data$Target) + 
                                        (max(data$Target) - min(data$Target)) * 0.2, 
                                      max(data$Target) + 
                                        (max(data$Target) - min(data$Target)) * 0.1)), 
               
               DiffRegPlus = sp + 
                 stat_cor(aes(color = Group), method = cor.method, size = 6, 
                          label.x = min(data$TF) + (max(data$TF) - min(data$TF)) * 0.05, 
                          label.y = c(max(data$Target) + 
                                        (max(data$Target) - min(data$Target)) * 0.2, 
                                      max(data$Target) + 
                                        (max(data$Target) - min(data$Target)) * 0.1)), 
               
               DysReg = sp + 
                 annotate("text",
                          label = paste("Regulatory intensity: ", 
                                        round(reg.i$unb.coef.1, 3), "; ", 
                                        paste(conf.int.level * 100, 
                                              "% CI: ", sep = ""),
                                        paste(round(reg.i$low.lim.1, 3), 
                                              round(reg.i$up.lim.1,3), 
                                              sep = "~"), 
                                        sep = ""), 
                          x = min(data$TF) + (max(data$TF) - min(data$TF)) * 0.4, 
                          y = max(data$Target) + (max(data$Target) - min(data$Target)) * 0.2, 
                          size = 5,colour = rainbow(5)[1]) + 
                 
                 annotate("text", 
                          label = paste("Regulatory intensity: ", 
                                        round(reg.i$unb.coef.2, 3), "; ", 
                                        paste(conf.int.level *  100, 
                                              "% CI: ", sep = ""), 
                                        paste(round(reg.i$low.lim.2, 3), 
                                              round(reg.i$up.lim.2, 3), 
                                              sep = "~"), 
                                        sep = ""), 
                          x = min(data$TF) + (max(data$TF) - min(data$TF)) * 0.4,
                          y = max(data$Target) + (max(data$Target) - min(data$Target)) * 0.1, 
                          size = 5, colour = rainbow(5)[4]))
  
  sp <- sp + theme_bw() + 
    theme(panel.background = element_rect(fill = "transparent"), 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          plot.background = element_rect(fill = "transparent"))
  
  sp <- sp + 
    theme(axis.title.x = element_text(face = "bold", size = 16), 
          axis.text.x = element_text(size = 16), 
          axis.title.y = element_text(face = "bold", size = 16), 
          axis.text.y = element_text(size = 16), 
          legend.title = element_text(size = 16,face = "bold"), 
          legend.text = element_text(size = 16))
  
  plot(sp)
}



################################################################################
#' @name DysregClin
#' @title Analyse the association of individual gene dysregualtion with clinical factor
#' @description Divide samples into different subgroups based the expression value of TF and target of a special gene dysregulation, and analyse its association with clinical factor.
#' 
#' @usage 
#' DysregClin(tf, tar, exp.data, 
#'            clin.data, clin.data.type, 
#'            divide.point = median)
#' 
#' @param tf The TF of a gene dysregulation.
#' @param tar The target of a gene dysregulation.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param clin.data The clinical data. The surv.data must include at least three colomns, "sample", "time", "status". 
#' @param clin.data.type The type of clinical data, such as "continuous", "discrete", or "survival".
#' @param divide.point The point used to divide samples into subgroups, such as median point.
#' 
#' @import survival
#' 
#' @return 
#' The significance of the association between gene dysregualtion and clinical factor. 
#' 
#' @examples 
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', 
#'                      pValue = 0.01, ci = 0.90, verbose = T)
#' 
#' dysreg.res <- dysreg.out$dysreg
#' 
#' 
#' # for continuous data type                                         
#' clin.data <- ClinData[,c("sample", "FMOne mutation burden per MB")]
#' clin.data <- clin.data[!is.na(clin.data$`FMOne mutation burden per MB`),]
#' head(clin.data)
#' 
#' DysregClin(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1], 
#'            exp.data = ExpData, clin.data = clin.data,
#'            clin.data.type = 'continuous')
#' 
#' 
#' #for discrete data type                                           
#' clin.data <- ClinData[,c("sample", "Immune phenotype")]
#' clin.data <- clin.data[!is.na(clin.data$`Immune phenotype`),]
#' head(clin.data)
#' 
#' DysregClin(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1], 
#'            exp.data = ExpData, clin.data = clin.data, 
#'            clin.data.type = 'discrete',
#'            divide.point = median)
#'                   
#' 
#' # for survival data type                                        
#' clin.data <- ClinData[,c("sample", "os","censOS")]
#' clin.data <- clin.data[!is.na(clin.data$os) & !is.na(clin.data$censOS),]
#' colnames(clin.data)[2:3] <- c('time','status')
#' head(clin.data)
#' 
#' DysregClin(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1], 
#'            exp.data = ExpData, clin.data = clin.data, 
#'            clin.data.type = 'survival')
#' }
#' 
#' @export
DysregClin <- function(tf, tar, exp.data, clin.data, clin.data.type, divide.point = median){
  
  if(clin.data.type == 'survival'){
    need.cols <- c("sample", "time", "status")
    judged.res <- is.element(need.cols, colnames(clin.data))
    if("FALSE" %in% judged.res) {
      stop("The input of clin.data must include at least three columns, that are 'sample', 'time', and 'status'. Please check the clin.data whether includes these columns and ensure the colname name is consistent with 'sample', 'times' and 'status'")
    }
  } else {
    need.cols <- "sample"
    judged.res <- is.element(need.cols, colnames(clin.data))
    if("FALSE" %in% judged.res){
      stop("The input of clin.data must include at least two columns, the first column is 'sample' the second column is clinical factor. Please check the clin.data whether includes the two columns")
    }
  }
  
  data <- data.frame(sample = colnames(exp.data),
                     TF = exp.data[rownames(exp.data) == tf,],
                     Tar = exp.data[rownames(exp.data) == tar,],
                     stringsAsFactors = F)
  
  divide.p <- apply(data[,c('TF','Tar')], 2, divide.point)
  
  data <- within(data,{
    group <- NA
    group[data$TF > divide.p[1] & data$Tar > divide.p[2]] <- 'TF_H & Tar_H'
    group[data$TF > divide.p[1] & data$Tar < divide.p[2]] <- 'TF_H & Tar_L'
    group[data$TF < divide.p[1] & data$Tar > divide.p[2]] <- 'TF_L & Tar_H'
    group[data$TF < divide.p[1] & data$Tar < divide.p[2]] <- 'TF_L & Tar_L'
  })
  
  if(clin.data.type == 'survival'){
    surv.data <- clin.data[,c('sample','time','status')]
    surv.data <- surv.data[!is.na(surv.data$time)& !is.na(surv.data$status),]
    surv.data$time <- as.numeric(surv.data$time)
    surv.data$status <- as.numeric(surv.data$status)
  
    data <- merge(data,surv.data,by = 'sample')
    
    ana.res <- survdiff(Surv(time, status) ~ group, data = data, rho = 0)
  }
  
  if(clin.data.type != 'survival'){
    
    colnames(clin.data)[2] <- 'clin.factor'
    clin.data <- clin.data[!is.na(clin.data$clin.factor),]
    
    data <- merge(data,clin.data,by = 'sample')
    
    if(clin.data.type == 'continuous'){
      sw.res <- shapiro.test(data$clin.factor)
      sw.res <- sw.res$p.value
      
      if(sw.res > 0.01){
        ana.res <- aov(clin.factor ~ group, data = data)
        ana.res <- list(method = 'ANOVA',res = ana.res)
      } else {
        ana.res <- kruskal.test(clin.factor ~ group, data = data)
        ana.res <- list(method = 'Kruskal-Wallis Rank Sum Test',res = ana.res)
      }
    }
    if(clin.data.type == 'discrete'){
      data <- table(data$group, data$clin.factor)
      ana.res <- chisq.test(data)
      ana.res <- list(method = 'Chi-squared Test',res = ana.res)

    }
  }
  
  return(ana.res)
}



################################################################################
#' @name plotDysregKM
#' @title Kaplan-Meier plots for a individual gene dysregualtion
#' @description Divide samples into different subgroups based the expression value of TF and target of a special gene dysregulation. plot Kaplan-Meier curves between the divided subgroups. 
#' 
#' @usage 
#' plotDysregKM(tf, tar, exp.data, surv.data, 
#'              divide.point = median, xlab, ylab, 
#'              pval = TRUE, pval.method = TRUE, 
#'              legend.position = c(0.8,0.9), ...)
#' 
#' @param tf The TF of a gene dysregulation.
#' @param tar The target of a gene dysregulation.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param surv.data The survival data. The surv.data  must include at least three colomns, "sample", "time", "status". 
#' @param divide.point The point used to divide samples into subgroups.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param pval Logical value. If TRUE, the p-value is added on the plot.
#' @param pval.method Whether to add a text with the test name used for calculating the pvalue corresponding to survival curves' comparison. Used only when pval=TRUE.
#' @param legend.position The position for presentating legend on the plot.
#' @param ... Other parameters passed to \link{ggsurvplot}.
#' 
#' @import survival
#' @import survminer
#' 
#' @examples 
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', 
#'                      pValue = 0.01, ci = 0.90, verbose = T)
#' 
#' surv.data <- ClinData[,c("sample", "os", "censOS")]
#' colnames(surv.data)[2:3] <- c('time','status')
#' 
#' plotDysregKM(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1], 
#'              exp.data = ExpData, surv.data = surv.data, 
#'              divede.point = median, 
#'              xlab = 'Months', ylab = 'Survival probability', 
#'              pval = TRUE, pval.method = TRUE,
#'              legend.position = c(0.8,0.9))
#'               
#' }
#' 
#' @export
plotDysregKM <- function(tf, tar, exp.data, surv.data, divide.point = median, 
                         xlab, ylab = 'Survival probability', 
                         pval = TRUE, pval.method = TRUE,
                         legend.position = c(0.8,0.9), ...){
    
  need.cols <- c("sample", "time", "status")
  judged.res <- is.element(need.cols, colnames(surv.data))
  if ("FALSE" %in% judged.res) {
    stop("The input of surv.data must include at least three columns, that are 'sample', 'time', and 'status'. Please check the surv.data whether includes these columns and ensure the colnames is consistent with 'sample', 'times' and 'status'")
  }
    
  data <- data.frame(sample = colnames(exp.data),
                     TF = exp.data[rownames(exp.data) == tf,],
                     Tar = exp.data[rownames(exp.data) == tar,],
                     stringsAsFactors = F)
  
  divide.p <- apply(data[,c('TF','Tar')],2,divide.point)
    
  data <- within(data,{
    group <- NA
    group[data$TF > divide.p[1] & data$Tar > divide.p[2]] <- 'TF_H & Tar_H'
    group[data$TF > divide.p[1] & data$Tar < divide.p[2]] <- 'TF_H & Tar_L'
    group[data$TF < divide.p[1] & data$Tar > divide.p[2]] <- 'TF_L & Tar_H'
    group[data$TF < divide.p[1] & data$Tar < divide.p[2]] <- 'TF_L & Tar_L'
  })
    
  surv.data <- surv.data[,c('sample','time','status')]
  surv.data <- surv.data[!is.na(surv.data$time)& !is.na(surv.data$status),]
  surv.data$time <- as.numeric(surv.data$time)
  surv.data$status <- as.numeric(surv.data$status)
    
  data <- merge(data,surv.data,by = 'sample')
  
  kmfit <- survfit(Surv(time, status) ~ group, data = data)
    
  p <- ggsurvplot(kmfit, data = data,
                  pval = TRUE,pval.method = TRUE,
                  legend.title = '',xlab = xlab, ylab=ylab,
                  font.main = c(16, "bold", "black"),
                  font.x = c(16,'bold','black'), 
                  font.y = c(16,'bold','black'),
                  font.tickslab = c(16, "plain", "black"))
    
  p <- p$plot + theme(legend.text = element_text(size = 16),
                      legend.position = legend.position)
    
  plot(p)
}


################################################################################
#' @name RankDysReg
#' @title Rank gene dysregulatins
#' @description Rank gene dysregualtions based on dysregulation degree quantified by combining differential regulation and differential expression.
#' @usage 
#' RankDysReg(dysreg)
#' 
#' @param dysreg The results of \code{\link{DysReg}}.
#' 
#' @return The rank of gene dysregulations.
#' 
#' @seealso 
#' \code{\link{RankDysTF}}
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', 
#'                      pValue = 0.01, ci = 0.90, verbose = T)
#'                      
#' dysreg.res <- dysreg.out$dysreg
#' 
#' reg.rank <- RankDysReg(dysreg = dysreg.res)
#' head(reg.rank)
#' 
#' }
#' 
#' @export
RankDysReg <- function(dysreg){
  
  net.i <- dysreg
  net.i$weight <- abs((net.i$unb.coef.2 - net.i$unb.coef.1) * net.i$de.logFC)
  dysreg.rank <- net.i[,c('TF','Target','weight')]
  
  dysreg.rank <- dysreg.rank[order(dysreg.rank$weight,decreasing = T),]
  
  rownames(dysreg.rank) <- NULL
  
  return(dysreg.rank)
}



################################################################################
#' @name RankDysTF
#' @title Rank TF with dysregulatin degree
#' @description Rank TF  based on dysregulation degree quantified by combining differential regulation and differential expression.
#' @usage 
#' RankDysTF(dysreg)
#' 
#' @param dysreg The results of \code{\link{DysReg}}.
#' 
#' @import igraph
#' 
#' @return The rank of TFs.
#' 
#' @seealso \code{\link{RankDysReg}}
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', 
#'                      pValue = 0.01, ci = 0.90, verbose = T)
#' 
#' dysreg.res <- dysreg.out$dysreg
#' 
#' tf.rank <- RankDysTF(dysreg = dysreg.res)
#' head(tf.rank)
#' 
#' }
#' 
#' @export

RankDysTF <- function(dysreg){
  
  net.i <- dysreg
  net.i$weight <- abs((net.i$unb.coef.2 - net.i$unb.coef.1) * net.i$de.logFC)
  net.i <- graph.data.frame(net.i[,c('TF','Target','weight')], directed = TRUE)
  
  degree <- strength(net.i) #weighted degree
  degree <- as.data.frame(degree)
  degree$Gene <- rownames(degree)
  tf.dys <- degree[degree$Gene %in% dysreg$TF,]
  
  tf.rank <- tf.dys[order(tf.dys$degree,decreasing = T),]
  
  tf.rank <- tf.rank[,c(2,1)]
  rownames(tf.rank) <- NULL
  
  return(tf.rank)
}



################################################################################
#' @name combineDysreg
#' @title Search the best combination of dysregulations for building signature
#' @description Combine dysregulations for building predictive signature with mechanistic interpretability by using genetic algorithm.
#' @usage 
#' combineDysreg(dysreg, exp.data,pheno.data, fitness.func, 
#'               pop.size = 1000, select.rate=0.2, add.rate=0.1, mut.rate = 0.1, 
#'               topN, train.rate = 0.6, iter=100, verbose = TRUE)
#' 
#' @param dysreg The dysregulations output from \code{\link{DysReg}}.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param fitness.func The  function for calculating fitness in genetic algorithm.
#' @param pheno.data The phenotype data corresponding to fitness.func.
#' @param pop.size The size of random populaiton generated in genetic algorithm.
#' @param select.rate The rate of selcting individual with good fitness from population.
#' @param add.rate The rate for adding individuals with not good fitness from population.
#' @param mut.rate The rate of mutation in genetic algorithm.
#' @param train.rate The rate of sample for training model in cross-validation.
#' @param topN The top N individuals with best fitness output in each iteration.
#' @param iter The times of iteration.
#' @param verbose A logical value indicating whether display the computing progress.
#' 
#' @seealso   
#' \code{\link{fitness.AUC}}; \code{\link{fitness.Cindex}}
#' 
#' @return 
#' The results of genetic algorithm:
#'  \item{individ}{The top N individuals with best fitness output in each iteration, from which useers could get the best combination of dysregulaitons and use it to build prdictive signatures.}
#'  \item{best.fitness}{The fitness value for top N individuals.}
#' 
#' @examples
#' \donttest{
#' data(ExpData)
#' data(tf2tar)
#' data(ClinData)
#' 
#' group.1 <- ClinData$sample[which(ClinData$binaryResponse == 'CR/PR')]
#' exp.1 <- ExpData[,colnames(ExpData) %in% group.1]
#' 
#' group.2 <- ClinData$sample[which(ClinData$binaryResponse == 'SD/PD')]
#' exp.2 <- ExpData[,colnames(ExpData) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(ExpData),1000)
#' 
#' tmp.1 <- exp.1[test.genes,]
#' tmp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = tmp.1, exp.2 = tmp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', 
#'                      pValue = 0.01, ci = 0.90, verbose = T)
#'                      
#' dysreg.res <- dysreg.out$dysreg
#' 
#' dysreg <- dysreg.res[,1:2]
#' 
#' pheno.data <- ClinData[,c("sample", "binaryResponse")]
#' pheno.data <- pheno.data[!is.na(pheno.data$binaryResponse),]
#' colnames(pheno.data)[2] <- "clin.factor"
#' head(pheno.data)
#' 
#' ## use genetic algorithm to search the best combination of dysregulations
#' combdysreg.out <- combineDysreg(dysreg = dysreg, exp.data = ExpData,
#'                                 fitness.func = 'fitness.AUC', 
#'                                 pheno.data = pheno.data, 
#'                                 pop.size = 1000, select.rate=0.2, 
#'                                 mut.rate=0.1, add.rate=0.1, 
#'                                 topN = 10, train.rate = 0.6, iter = 100)
#' 
#' ## Check the output of combineDysreg
#' best.fit <- combdysreg.out$best.fitness
#' best.individ <- combdysreg.out$individ
#' head(best.fit)
#' summary(best.fit$fitness)
#' 
#' best.individ[which.max(best.fit$fitness),]
#' 
#' ## THe example for building signature from output of combineDysreg
#' 
#' # calculate the frequence of each each dysregualtion among the top N best individuals
#' for(i in 1:100){
#'   start <- 10*i-9
#'   end <- 10*i
#'   iter.i <- best.individ[start:end,]
#'   iter.i <- colSums(iter.i)/10
#'   freq.res <- rbind(freq.res,iter.i)
#' }
#' 
#' freq.res <- t(freq.res)
#' colnames(freq.res) <- c(1:100)
#' rownames(freq.res) <- paste(dysreg$TF,dysreg$Target,sep = '-')
#' freq.res <- freq.res[order(rownames(freq.res)),]
#' 
#' # visualize the frequence of each each dysregualtion among the top N best individuals
#' pheatmap(freq.res,color = colorRampPalette(c('lightyellow',"orange", "firebrick3"))(1000), 
#'          display_numbers = F, cluster_rows = F,cluster_cols = F,
#'          fontsize_row = 8, fontsize_col = 10)
#'          
#' # choose the dysregulations frequently emerged among the top N best individuals as signatures
#' markers <- rownames(freq.res)[freq.res[,100] >= 0.9]
#' mark.genes <- unique(unlist(strsplit(markers,'-')))
#' 
#' }
#' 
#' @export
combineDysreg <- function(dysreg, exp.data, pheno.data, fitness.func, 
                          pop.size = 1000, select.rate = 0.2,add.rate = 0.1, 
						              mut.rate = 0.1, topN = 10,
                          train.rate = 0.6,iter = 100, verbose = TRUE){
  
  if (fitness.func == 'fitness.Cindex'){
    need.cols <- c("sample", "time", "status")
    judged.res <- is.element(need.cols, colnames(pheno.data))
    
    if ("FALSE" %in% judged.res) {
      stop("The input of pheno.data must include at least three columns, that are 'sample', 'time', and 'status'. Please check the pheno.data whether includes these columns and ensure the colnames are consistent with 'sample', 'times' and 'status'")
    }
    
    pheno.data$time <- as.numeric(pheno.data$time)
    pheno.data$status <- as.numeric(pheno.data$status)
    pheno.data <-pheno.data[!is.na(pheno.data$time) & !is.na(pheno.data$status),]
  } else {
    need.cols <- "sample"
    judged.res <- is.element(need.cols, colnames(pheno.data))
    if ("FALSE" %in% judged.res) {
      stop("The input of pheno.data must include at least two columns, one column is 'sample' the other column is clinical factor. Please check the pheno.data whether includes the two columns")
    }
    colnames(pheno.data)[2] <- 'clin.factor'
    pheno.data <- pheno.data[!is.na(pheno.data$clin.factor),]
  }
  
  chr.len <- nrow(dysreg)
  if(chr.len > nrow(pheno.data)*0.5){
    stop('The number of inputting dysregulations is larger than half of example size, please filter dysregulations.')
  }
  
  ## generate a random initial population
  pop <- do.call(rbind,lapply(1:pop.size,function(i,chr.len){
    gene.num <- sample(1:chr.len,1,replace = F)
    individ <- sample(c(rep(1,gene.num),rep(0, chr.len - gene.num)),
                      chr.len,replace = F)
    return(individ)
  },chr.len = chr.len))
  
  
  bestfit.out <- list()
  bestfit.out$individ <- vector()
  bestfit.out$best.fitness <- vector()
  
  for (i in 1: iter) {
    
    rownames(pop) <- c(1:nrow(pop))
    
    ##calculate fitness
    obj.fitness <- suppressMessages(apply(pop, 1, fitness.func,
                                    dysreg = dysreg,exp.data = exp.data,
                                    pheno.data = pheno.data,
                                    train.rate = train.rate))
    
    obj.fitness <- as.data.frame(obj.fitness)
    obj.fitness$indival <- rownames(obj.fitness)
    obj.fitness <- obj.fitness[,c(2,1)]
    colnames(obj.fitness)[2] <- 'obj.fitness'
    
    pop.obj <- cbind(pop,obj.fitness)
    pop.obj <- pop.obj[order(pop.obj$obj.fitness,decreasing = T),]
    obj.value <- obj.fitness[order(obj.fitness$obj.fitness,decreasing = T),]
    
    best.fitness <- data.frame(iter = i,
                               fitness = obj.value$obj.fitness[1:topN],
                               stringsAsFactors = F)
    
    bestfit.out$best.fitness <- rbind(bestfit.out$best.fitness,
                                      best.fitness)
    
    best.individ <- pop.obj[pop.obj$indival %in% obj.value$indival[1:topN],
                            1:(ncol(pop.obj)-2)]
    bestfit.out$individ <- rbind(bestfit.out$individ,best.individ)
    
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = iter, style = 3)
      setTxtProgressBar(pb,i)
    }
    
    ## evolution prorecess
    if(i < (iter-1)){
      children.pop <- evolution(pop = pop, pop.obj = pop.obj, 
                                obj.value = obj.value, 
                                pheno.data = pheno.data,
                                pop.size = pop.size, select.rate = 0.2, 
                                mut.rate = 0.1,add.rate = 0.1)
      
      pop <- children.pop
    }
    
  }
  
  return(bestfit.out)
}


evolution <- function(pop, pop.obj, obj.value, 
                      pop.size = pop.size, pheno.data,
                      select.rate = select.rate, add.rate = add.rate, 
                      mut.rate = mut.rate){
  
  ##select
  selected.individ <- obj.value[1:round((nrow(pop) * select.rate)),1]
  deleted.individ <- obj.value[-(1:round((nrow(pop) * select.rate))),1]
  
  added.individ <- deleted.individ[runif(length(deleted.individ)) < add.rate]
  selected.individ <- c(selected.individ,added.individ)
  
  selected.pop <- pop.obj[pop.obj$indival %in% selected.individ,]
  selected.obj <- obj.value[obj.value$indival %in% selected.individ,]
  
  ##crossover
  desired.len <- pop.size
  
  children.pop <- vector()
  
  for (i in 1:(nrow(selected.pop)-1)){
    individ.1 <- as.numeric(selected.pop[i,1:(ncol(selected.pop)-2)])
    individ.2 <- as.numeric(selected.pop[i+1,1:(ncol(selected.pop)-2)])
    
    cross.point <- round(length(individ.1)*runif(1,min = 0.1,max = 0.9))
    
    child.1 <- c(individ.1[1:cross.point],individ.2[-(1:cross.point)])
    child.2 <- c(individ.2[1:cross.point],individ.1[-(1:cross.point)])
    child2 <- rbind(child.1,child.2)
    
    ##mutation
    random.rate <- runif(nrow(child2))
    if(length(random.rate[random.rate < mut.rate])>0){
      mut.individ <- which(random.rate<mut.rate)
      for(m in mut.individ){
        pos.mutate <- sample(1:ncol(child2),1)
        child2[m,pos.mutate] <- abs(child2[m,pos.mutate] - 1)
      }
    }
    
    gene.num <- apply(child2,1,sum)
    child2 <- child2[gene.num > 0 & gene.num <= nrow(pheno.data)*0.5,]
    
    children.pop <- rbind(children.pop,child2)
  }
  
  children.pop.len <- nrow(children.pop) 
  
  while(children.pop.len < desired.len){
    selected.obj$prob <- exp(selected.obj$obj.fit)/sum(exp(selected.obj$obj.fit))
    individ.1 <- sample(selected.obj$indival[-1],1,prob = selected.obj$prob[-1])
    individ.2 <- selected.obj$indival[which(selected.obj$indival==individ.1)-1]
    
    individ.1 <- as.numeric(selected.pop[selected.pop$indival==individ.1,
                                         1:(ncol(selected.pop)-2)])
    individ.2 <- as.numeric(selected.pop[selected.pop$indival==individ.2,
                                         1:(ncol(selected.pop)-2)])
    
    cross.point <- round(length(individ.1)*runif(1,min = 0.1,max = 0.9))
    child.1 <- c(individ.1[1:cross.point],individ.2[-(1:cross.point)])
    child.2 <- c(individ.2[1:cross.point],individ.1[-(1:cross.point)])
    child2<-rbind(child.1,child.2) 
    
    ##mutation
    random.rate <- runif(nrow(child2))
    if(length(random.rate[random.rate < mut.rate])>0){
      mut.individ <- which(random.rate<mut.rate)
      for(m in mut.individ){
        pos.mutate <- sample(1:ncol(child2),1)
        child2[m,pos.mutate] <- abs(child2[m,pos.mutate] - 1)
      }
    }
    
    ## new generateion
    gene.num <- apply(child2,1,sum)
    child2 <- child2[gene.num > 0 & gene.num <= nrow(pheno.data)*0.5,]
    
    children.pop <- rbind(children.pop,child2)
    
    children.pop.len <- nrow(children.pop)  
  }
  
  return(children.pop)
}



################################################################################
#' @name fitness.Cindex
#' @title Calculate C-Index within survival data for using as fitness function
#' @description Calculate C-Index with cross-validation in survival data for using as fitness function.
#' @usage 
#' fitness.Cindex(combine.dysreg, dysreg, exp.data, pheno.data, train.rate)
#' 
#' @param combine.dysreg The combination of dysregulations.
#' @param dysreg The dysregulations output from \code{\link{DysReg}}.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param pheno.data The phenotype data corresponding to fitness function. At here, it is survival data.
#' @param train.rate The rate of sample for training model in cross-validation.
#' 
#' @import survival
#' @import survcomp
#' 
#' @return 
#' The fitness, C-Index, for each individual.
#' 
#' @seealso 
#' \code{\link{combineDysreg}}
#' 
#' @export
fitness.Cindex <- function(combine.dysreg,dysreg,exp.data,pheno.data,train.rate){
	
  regs <- dysreg[combine.dysreg !=0,]
  genes <- union(regs[,1],regs[,2])
  
  exp <- exp.data[rownames(exp.data) %in% genes,]
  exp <- as.data.frame(t(exp),stringsAsFactors = F)
  exp$sample <- rownames(exp)
  colnames(exp) <- sub('-','_',colnames(exp))
  
  colnames(pheno.data)[2:3] <- c('time','status')
  pheno.data <- pheno.data[!is.na(pheno.data$time) & !is.na(pheno.data$status),] 
  
  data <- merge(exp, pheno.data, by = 'sample')
  data <- data[,-1]
  
  cv.ci <- do.call(rbind,lapply(1:10,function(i,data, train.rate){
    cv.select <- sample(1:nrow(data),round(nrow(data)*train.rate),replace = F)
    train.data <- data[cv.select,]
    test.data <- data[setdiff(1:nrow(data),cv.select),]
    
    genes <- setdiff(colnames(exp),'sample')
    ff <- as.formula(paste("Surv(time,status)~", paste(genes, collapse = "+")))
    train.mod <- coxph(ff, data=train.data)
    train.mod$coefficients[is.na(train.mod$coefficients)] <- 0 
    score <- as.matrix(test.data[,1:(ncol(test.data) - 2)])%*% train.mod$coefficients
    res.i <- data.frame(iter = i,
                        fitness = concordance.index(score,test.data$time,test.data$status)$c.index)
    return(res.i)
  }, data = data, train.rate = train.rate))
  
  cv.ci <- mean(cv.ci$fitness,trim = 0.1)
  return(cv.ci)
}



################################################################################
#' @name fitness.AUC
#' @title Calculate AUC in classifying data for using as fitness function
#' @description Calculate AUC with cross-validation in classifying data for using as fitness function.
#' @usage 
#' fitness.Cindex(combine.dysreg, dysreg, exp.data, pheno.data, train.rate)
#' 
#' @param combine.dysreg The combination of genes.
#' @param dysreg The dysregulations output from \code{\link{DysReg}}.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param pheno.data The phenotype data corresponding to objective function. At here, it is classifying data.
#' @param train.rate The rate of sample for training model in cross-validation.
#' 
#' @import ROCR
#' @import pROC
#' @import e1071
#' 
#' @return 
#' The fitness, AUC, for each individual.
#' 
#' @seealso 
#' \code{\link{combineDysreg}}
#' 
#' @export
fitness.AUC <- function(combine.dysreg, dysreg, exp.data, pheno.data, train.rate) {
  regs <- dysreg[combine.dysreg != 0, ]
  genes <- union(regs[, 1], regs[, 2])
  exp <- exp.data[rownames(exp.data) %in% genes, ]
  exp <- as.data.frame(t(exp), stringsAsFactors = F)
  exp$sample <- rownames(exp)
  colnames(exp) <- sub("-", "_", colnames(exp))
  
  pheno.data <- pheno.data[!is.na(pheno.data$clin.factor), ]
  
  data <- merge(exp, pheno.data, by = "sample")
  data <- data[, -1]
  
  cv.auc <- do.call(rbind, lapply(1:10, function(i, data, train.rate){
    cv.select <- sample(1:nrow(data), round(nrow(data) * 
                                              train.rate), replace = F)
    train.data <- data[cv.select, ]
    test.data <- data[setdiff(1:nrow(data), cv.select), ]
    genes <- setdiff(colnames(exp), c('clin.factor'))
    
    x <- train.data[,colnames(train.data) %in% genes]
    y <- as.factor(train.data$clin.factor)
    model <- svm(x, y)
    pred <- predict(model, test.data[,colnames(test.data) %in% genes], 
                    decision.values = TRUE)
    test.data$score=attr(pred, "decision.values")
    
    predob<- prediction(test.data$score, test.data$clin.factor)
    perf<- performance(predob, 'tpr','fpr')
    auc.ci=signif(ci(test.data$clin.factor,test.data$score)[1:3],3)
    
    res.i <- data.frame(times = i, fitness=auc.ci[2],stringsAsFactors = F)
    
    return(res.i)
  }, data = data, train.rate = train.rate))
  
  cv.auc <- mean(cv.auc$fitness, trim = 0.1)
  return(cv.auc)
}



################################################################################
#' @name ExpData
#' @docType data
#' @title Example of expression data
#' @description An example of expression data.
#' @usage data(ExpData)
#' 
#' @details The example data was derived from IMvigor210CoreBiologies, which includes RNA-seq data of 348 samples and corresponding drug response of atezolizumab.
#' 
#' @references 
#' Mariathasan S, Turley S J, Nickles D, et al. TGFbeta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018, 554(7693): 544-8.
#'
#' @examples
#' data(ExpData)
#' dim(ExpData)
#' ExpData[1:5,1:5,]
#' 
#' @keywords datasets
#' @export



################################################################################
#' @name ClinData
#' @docType data
#' @title Example of clinical data
#' @description An example of clinical data.
#' @usage data(ClinData)
#' 
#' @details The example data was derived from IMvigor210CoreBiologies, which includes RNA-seq data of 348 samples and corresponding drug response of atezolizumab.
#' 
#' @references 
#' Mariathasan S, Turley S J, Nickles D, et al. TGFbeta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018, 554(7693): 544-8.
#' 
#' @examples
#' data(ClinData)
#' head(ClinData)
#' 
#' @keywords datasets
#' @export


################################################################################
#' @name tf2tar
#' @docType data
#' @title TF-target regulations
#' @description TF-target regulations within a reference GRN.
#' @usage data(tf2tar)
#' 
#' @details The TF-target relationships were derived from predicting by FIMO with TF motif data, which was downloaded from HumanTF database.
#' 
#' @references 
#' Grant C E, Bailey T L, Noble W S. FIMO: scanning for occurrences of a given motif. Bioinformatics. 2011, 27(7): 1017-8.
#' @references 
#' Lambert S A, Jolma A, Campitelli L F, et al. The Human Transcription Factors. Cell. 2018, 172(4): 650-65.
#' 
#' @examples
#' data(tf2tar)
#' head(tf2tar)
#' 
#' @keywords datasets
#' @export


################################################################################
#' @name DE
#' @docType data
#' @title Example of DEGs input for function DysReg and DiffRegPlus
#' @description Example of DEGs input for function DysReg and DiffRegPlus.
#' @usage data(DE)
#' 
#' @details de.genes contains three columns, GeneSymbol is the gene symbol of DEGs, high.condition indicates which comdition expresses high level, de.logFC is logFC value.
#' 
#' @examples
#' data(DE)
#' head(DE)
#' 
#' @keywords datasets
#' @export