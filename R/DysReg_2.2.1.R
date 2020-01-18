
####################################################################################################
#' @name condiGRN
#' @title Build conditional GRN
#' @description Build conditional gene regulatory network (GRN) with expression data and prior network by using feature selection algorithm.
#' @usage condiGRN(exp, tf2tar, method , pValue, threshold)
#' 
#' @param exp Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param tf2tar The TF-target interaction pairs within a prior GRN.
#' @param method The method used,such as 'Boruta', 'RGBM'.
#' @param pValue Confidence level used in Boruta.
#' @param threshold The threshould for weight in RGBM.
#' 
#' @import Boruta
#' @import RGBM
#' @import igraph
#' 
#' @details While using Boruta, the threshold is the predifined pValue. While using RGBM or GENIE, users could set threshould based on output of weight.
#' 
#' @return A conditional GRN.
#' 
#' @examples
#' \donttest{
#' # Build a conditional GRN based on a prior GRN.
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' data(clin_data)
#' head(clin_data)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' temp <- exp.1[sample(1:nrow(exp.1),1000),]
#' 
#' ## using method Boruta
#' net.1 <- condiGRN(exp = temp, tf2tar = tf2tar, method = 'Boruta', pValue = 0.01)
#' 
#' }
#' 
#' @export
condiGRN <- function(exp, tf2tar, method = 'Boruta', pValue = 0.01, threshold){
  
  if (ncol(exp) < 5) {
    stop("Expression data must have at least five samples.")
  }
  if (ncol(exp) >= 5 & ncol(exp) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp[is.na(exp)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  tf2tar <- tf2tar[tf2tar$Target %in% rownames(exp),]
  tf2tar <- tf2tar[tf2tar$TF %in% rownames(exp),]
  
  ## build conditinal GRN with Boruta
  if (method == 'Boruta'){
    targets <- unique(tf2tar$Target)
    
    condition.GRN <- do.call(rbind,lapply(1:length(targets),function(i,exp,tf2tar){
      reg.tar.i <- tf2tar[tf2tar$Target == targets[i],]
      reg.i <- unique(reg.tar.i$TF)
      
      if(length(reg.i) > 1){
        exp.i <- as.data.frame(t(exp[c(targets[i],reg.i),]))
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
      
      
      return(reg.tar.i)
      
    },exp = exp,tf2tar = tf2tar))
  }
  
  ## build conditinal GRN with RGBM
  if (method == 'RGBM'){
   
    exp <- t(exp)
    
    net.matrix <- graph.data.frame(tf2tar)
    net.matrix <- as_adjacency_matrix(net.matrix)
    
    K <- matrix(0, nrow(exp), ncol(exp))
    colnames(K) <- colnames(exp)
    rownames(K) <- rownames(exp)
    
    condition.GRN <- RGBM(E = exp,K = K,
                          g_M = net.matrix, 
                          tfs = unique(tf2tar$TF),
                          targets = unique(tf2tar$Target))
    
    condition.GRN <- melt(condition.GRN, na.rm = TRUE)
    colnames(condition.GRN) <- c("TF", "Target", "weight")
    condition.GRN <- condition.GRN[condition.GRN$weight > threshold, ]
  }
  
  condition.GRN$TF <- as.character(condition.GRN$TF)
  condition.GRN$Target <- as.character(condition.GRN$Target)
  
  rownames(condition.GRN) <- NULL
  
  return(condition.GRN)
} 



####################################################################################################
#' @name quantiReg
#' @title Quantify regulatory intensities of regulations
#' @description Quantify regulatory intensities with de-biased LASSO.
#' @usage quantiReg(exp, net, ci)
#' 
#' @param exp Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param net Conditioanl GRN outpou from \code{condiGRN}.
#' @param ci  Confident invetals of coefficients.
#'
#' @import 
#'
#' @return A data frame containing the regulatory intenty and its confident invertal for each regulation.
#' 
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' data(clin_data)
#' head(clin_data)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' temp <- exp.1[sample(1:nrow(exp.1),100),]
#' 
#' ## Build a conditional GRN based on a prior GRN
#' net.1 <- condiGRN(exp = temp, tf2tar = tf2tar, method = 'Boruta', pValue = 0.01)
#' 
#' ## Quantify regulatory intensity
#' quanti.net.1 <- quantiReg(exp = temp, net = net.1, ci = 0.95)
#' }
#' 
#' @export
quantiReg <- function(exp, net, ci = 0.95){
  
  #source('lasso_inference.r')
  
  if (ncol(exp) < 5) {
    stop("Expression data must have at least five samples.")
  }
  if (ncol(exp) >= 5 & ncol(exp) < 10) {
    warning("Expression data have less than ten samples.")
  }
  
  if(length(exp[is.na(exp)]) > 0){
    stop("Please deal with the missing value in expression data.")
  }
  
  
  ## quantify regulatory intensity with de-baised LASSO
  net <- net[net$Target %in% rownames(exp),]
  net <- net[net$TF %in% rownames(exp),]
  
  targets <- unique(net$Target)
  
  quanti.reg <- do.call(rbind,lapply(1:length(targets),function(i,exp,net,ci){
    
    reg.tar.i <- net[net$Target == targets[i],]
    reg.i <- unique(reg.tar.i$TF)
    
    exp.i <- as.data.frame(t(exp[c(targets[i],reg.i),]),stringsAsFactors = F)
    colnames(exp.i)[1] <- 'tar'
    
    exp.i <- scale(exp.i,center =TRUE,scale = TRUE)
    
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
      glm.fit.ci <- confint(glm.fit, level = ci,method = 'confint.glm')
      
      ci.i <- data.frame(TF = rownames(glm.fit.coef)[-1], Target = targets[i],
                         unb.coef = glm.fit.coef[-1,1],
                         low.lim = glm.fit.ci[-1,1], up.lim = glm.fit.ci[-1,2],
                         P.val = glm.fit.coef[-1,4], stringsAsFactors = F)
      
    }
    
    return(ci.i)
    
  },exp = exp,net = net,ci = ci))
  
  rownames(quanti.reg)=NULL
  
  return(quanti.reg)
  
}



####################################################################################################
#' @name DysReg
#' @title Identify gene dysregulations
#' @description Identify gene dysregulations by integrating three properties including differential regulation, differential expression of target, and the consistency between differential regulation and differential expression.
#' @usage DysReg(exp.1, exp.2, tf2tar, de.genes = NULL,
#'               de.pval = NULL, de.qval = NULL, de.logFC = NULL, 
#'               grn.method = 'Boruta', 
#'               pValue = 0.01, threshold = NULL, ci = 0.95)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param de.genes A dataframe for differential expression genes. The dataframe must include three colomns, "GeneSymbol", "high.condition", "logFC". "high.condition" means which condition represents high expression. "logFC" is the output logFC from differential expression analysis.
#' @param tf2tar The TF-target interaction pairs within a prior GRN.
#' @param de.pval The cutoff of pval used for filtering differential expression genes.
#' @param de.qval The cutoff of qval used for filtering differential expression genes.
#' @param de.logFC The cutoff of absolute logFC used for filtering differential expression genes.
#' @param pValue Confidence level used in Boruta.
#' @param grn.method The method used for constructing GRN,such as 'Boruta', 'RGBM'. 
#' @param threshold The threshould for weight in RGBM.
#' @param ci  The confident invetals of coefficients.
#' 
#' @details DysReg first build conditional GRNs with a random forest-based feature selection algorithm, where each link??s regulatory intensity and its confidential interval is estimated with a de-biased LASSO method. Gene dysregulations were then identified by integrating three properties including differential regulation, differential expression of target, and the consistency between differential regulation and differential expression.
#' @import limma
#' @import Boruta
#' @import RGBM
#' 
#' @return 
#' The results of gene dysregulation analysis:
#'  \item{de.genes}{The identified differential expression genes between conditions.}
#'  \item{net.1}{The conditional GRN for exp.1.}
#'  \item{quanti.reg.1}{The quantification of regulatory intensity for exp.1.}
#'  \item{net.2}{The conditional GRN for exp.2.}
#'  \item{quanti.reg.2}{The quantification of regulatory intensity for exp.2.}
#'  \item{dysreg}{The identified gene dysregulations}
#'  
#'  @seealso condiGRN, quantiReg.
#'  
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' data(clin_data)
#' head(clin_data)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[,colnames(exp) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(exp),1000)
#' 
#' temp.1 <- exp.1[test.genes,]
#' temp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = temp.1, exp.2 = temp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', pValue = 0.01, ci = 0.90)
#' }
#' 
#' @export
DysReg <- function(exp.1, exp.2, tf2tar, de.genes = NULL,
                   de.pval = NULL, de.qval = NULL, de.logFC = NULL, 
                   grn.method = 'Boruta', pValue = 0.01, threshold = NULL, ci = 0.95){
  
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
    warning("Please check the high expression condition in your datasets.")
    }
    
  if(is.null(de.genes)){
    cat('\n',paste('Identifying differential expression genes'),'\n')
    
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
    contrast.matrix <- makeContrasts(paste0(c('condition1','condition2'),collapse = "-"),levels = design)
    
    fit <- lmFit(merged.exp,design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    diff.out <- topTable(fit2, coef = 1, n = Inf)
    diff.out <- na.omit(diff.out)
    
    diff.filter <- diff.out
    colnames(diff.filter)[c(1,4,5)] <- c('de.logFC','de.pval','de.qval')
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
      
      de.filter <- switch (cutoff.type,
                           param_0 = diff.filter,
                           param_1 = subset(diff.filter, diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,2]),
                           param_2 = subset(diff.filter, diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,2] & 
                                              diff.filter[,match(cutoff.s[2],colnames(diff.filter))] < cutoff[,3]))
    } else {
      cutoff.type <- paste('param', length(cutoff.s), sep = '_')
      
      de.filter <- switch (cutoff.type,
                           param_0 = diff.filter,
                           param_1 = subset(diff.filter, diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,1]),
                           param_2 = subset(diff.filter, diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,1] & 
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
    cat('\n',paste('The number of DEGs:',nrow(de.genes)),'\n')
  }
  
  cat('\n',paste('Modelling conditional GRN for exp.1'),'\n')
  
  tf2tar.1 <- tf2tar[tf2tar$TF %in% rownames(exp.1),]
  tf2tar.1 <- tf2tar.1[tf2tar.1$Target %in% rownames(exp.1),]
  tf2tar.1 <- tf2tar.1[tf2tar.1$Target %in% de.genes$GeneSymbol,]
  
  net.1 <- condiGRN(exp = exp.1,tf2tar = tf2tar.1, method = grn.method, 
                    pValue = pValue, threshold = threshold)
  
  
  cat('\n',paste('Modelling conditional GRN for exp.2'),'\n')
  
  tf2tar.2 <- tf2tar[tf2tar$TF %in% rownames(exp.2),]
  tf2tar.2 <- tf2tar.2[tf2tar.2$Target %in% rownames(exp.2),]
  tf2tar.2 <- tf2tar.2[tf2tar.2$Target %in% de.genes$GeneSymbol,]
  
  net.2 <- condiGRN(exp = exp.2,tf2tar = tf2tar.2, method = grn.method, 
                    pValue = pValue, threshold = threshold)
  
  net <- unique(rbind(net.1,net.2))
  
 
  
  cat('\n',paste('Quantifying regulatory intensity for exp.1'), '\n')
  quanti.reg.1 <- quantiReg(exp = exp.1,net = net,ci = ci)
  
  cat('\n',paste('Quantifying regulatory intensity for exp.2'), '\n')
  quanti.reg.2 <- quantiReg(exp = exp.2,net = net,ci = ci)
  
  
  dysreg.net <- merge(quanti.reg.1, quanti.reg.2, by = c('TF','Target'),all = T)
  dysreg.net[,3:10][is.na(dysreg.net[,3:10])] <- 0
  
  data.1 <- dysreg.net[dysreg.net$unb.coef.x * dysreg.net$unb.coef.y != 0,]
  data.2 <- dysreg.net[dysreg.net$unb.coef.x * dysreg.net$unb.coef.y == 0,]
  
  data.1 <- data.1[data.1$low.lim.x > data.1$up.lim.y | data.1$up.lim.x < data.1$low.lim.y,]
  data.2 <- data.2[data.2$low.lim.x * data.2$up.lim.x > 0 | data.2$low.lim.y * data.2$up.lim.y > 0,]
  
  dysreg <- rbind(data.1,data.2)
  rm(data.1)
  rm(data.2)
  
  dysreg <- merge(dysreg,de.genes,by.x = 'Target',by.y = 'GeneSymbol')
  dysreg <- dysreg[,c(2,1,3:12)]
  
  dysreg <- rbind(subset(dysreg, dysreg$high.condition == 1 & (dysreg$unb.coef.y - dysreg$unb.coef.x) < 0),
                  subset(dysreg, dysreg$high.condition == 2 & (dysreg$unb.coef.y - dysreg$unb.coef.x) > 0))
  
  dysreg.res <- list(de.genes = de.genes, net.1 = net.1, quanti.reg.1 = quanti.reg.1,
                     net.2 = net.2, quanti.reg.2 = quanti.reg.2,
                     dysreg = dysreg)
  
  return(dysreg.res)
}




####################################################################################################
#' @name DiffCor
#' @title Differential correlation analysis
#' @description Differential correlation analysis with Fishers' Z test.
#' @usage DiffCor(exp.1, exp.2, tf2tar, cor.method = 'pearson', p.adj)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param tf2tar The TF-target interaction pairs within a prior GRN.
#' @param cor.method Which correlation coefficient (or covariance) is to be computed. One of "pearson" (default) or "spearman", can be abbreviated.
#' @param p.adj  Correction method for p value adjust.
#' 
#' @return The identified differential correlation regulations.
#' 
#' @examples
#' \donttest{
#' 
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[,colnames(exp) %in% group.2]
#' 
#' dc.res <- DiffCor(exp.1, exp.2, tf2tar, 
#'                   cor.method = 'pearson', p.adj = 'BH')
#'                    
#' dc.res <- dc.res[dc.res$p.adj < 0.05,]
#' head(diff.res)
#' }
#' 
#' @export
DiffCor <- function(exp.1, exp.2, tf2tar, 
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
  
  tf2tar <- tf2tar[tf2tar$TF %in% merged.genes,]
  tf2tar <- tf2tar[tf2tar$Target %in% merged.genes,]
  
  tf.list=unique(tf2tar$TF)
  
  cor.1=cor(t(exp.1),method = cor.method)
  cor.1=cor.1[colnames(cor.1) %in% tf.list,]
  
  cor.2=cor(t(exp.2),method = cor.method)
  cor.2=cor.2[colnames(cor.2) %in% tf.list,]
  
  z.1=(0.5*log((1 + cor.1)/(1 - cor.1)))
  z.2=(0.5*log((1 + cor.2)/(1 - cor.2)))
  
  z.value=(z.2 - z.1)/((1/(ncol(exp.1)-3) + 1/(ncol(exp.2)-3))^0.5)
  
  diffcor.tf2tar=do.call(rbind,lapply(1:length(tf.list), function(i, tf2tar, tf.list, cor.1, cor.2, z.value){
    
    cor1.i <- cor.1[rownames(cor.1) == tf.list[i],]
    cor1.i <- as.data.frame(cor1.i)
    cor1.i$TF=tf.list[i]
    cor1.i$Target=rownames(cor1.i)
    
    cor2.i <- cor.2[rownames(cor.2) == tf.list[i],]
    cor2.i <- as.data.frame(cor2.i)
    cor2.i$TF=tf.list[i]
    cor2.i$Target=rownames(cor2.i)
    
    diffcor.i=merge(cor1.i,cor2.i,by = c('TF','Target'))
    
    z.i=z.value[rownames(z.value)==tf.list[i],]
    z.i=as.data.frame(z.i)
    z.i$TF=tf.list[i]
    z.i$Target=rownames(z.i)
    
    diffcor.i=merge(diffcor.i,z.i,by = c('TF','Target'))
    
    
    tf.i=tf2tar[tf2tar$TF == tf.list[i],]
    diffcor.i=diffcor.i[diffcor.i$Target %in% tf.i$Target,]
    colnames(diffcor.i)[3:5]=c('cor.1','cor.2','z.value')
    rownames(diffcor.i)=NULL
    
    return(diffcor.i)
  }, tf2tar=tf2tar,tf.list=tf.list,cor.1 = cor.1, cor.2=cor.2, z.value=z.value))
  
  diffcor.tf2tar <- within(diffcor.tf2tar, {
    p.val=0
    p.val[z.value > 0]= pnorm(z.value[z.value > 0],lower.tail = F)
    p.val[z.value < 0]= pnorm(z.value[z.value < 0],lower.tail = T)
  })
  
  diffcor.tf2tar$p.adj=p.adjust(diffcor.tf2tar$p.val,method = p.adj)
  
  return(diffcor.tf2tar) 
}



####################################################################################################
#' @name DiffCorPlus
#' @title Differential correlation analysis plus differential expression of target, and the consistency between differential correlation and differential expression.
#' 
#' @description Identify regulations with differential correlation analysis, differential expression of target, and the consistency between differential correlation and differential expression.
#' 
#' @usage DiffCorPlus(exp.1, exp.2, tf2tar, de.genes, cor.method = 'pearson', p.adj,
#'                    de.pval=NULL, de.qval = NULL, de.logFC = NULL)
#' 
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param tf2tar The TF-target interaction pairs within a prior GRN.
#' @param de.genes A dataframe for differential expression genes. The dataframe must include three colomns, "GeneSymbol", "high.condition", "logFC". "high.condition" means which condition represents high expression. "logFC" is the output logFC from differential expression analysis.
#' @param cor.method Which correlation coefficient (or covariance) is to be computed. One of "pearson" (default) or "spearman", can be abbreviated.
#' @param p.adj  Correction method for p value adjust.
#' @param de.pval The cutoff of pval used for filtering differential expression genes.
#' @param de.qval The cutoff of qval used for filtering differential expression genes.
#' @param de.logFC The cutoff of absolute logFC used for filtering differential expression genes.
#' 
#' @import limma
#' 
#' @return The identified regulations.
#' 
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[,colnames(exp) %in% group.2]
#' 
#' diffcorpp.res <- DiffCorPlus(exp.1,exp.2, tf2tar, de.genes = NULL, 
#'                              cor.method = 'pearson', p.adj = 'BH',
#'                              de.pval=0.05)
#' 
#' diffcorpp.res <- diffcorpp.res[diffcorpp.res$p.val < 0.05,]
#' head(diffcorpp.res)
#' }
#' 
#' @export
DiffCorPlus <- function(exp.1,exp.2, tf2tar,  de.genes = NULL, 
                         cor.method = 'pearson', p.adj = 'BH',
                         de.pval=NULL, de.qval = NULL, de.logFC = NULL){
  
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
    warning("Please check the high expression condition in your datasets.")
  }
  
  if(is.null(de.genes)){
    
    cat('\n',paste('Identifying differential expression genes'),'\n')
    
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
    contrast.matrix <- makeContrasts(paste0(c('condition1','condition2'),collapse = "-"),levels = design)
    
    fit <- lmFit(merged.exp,design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    diff.out <- topTable(fit2, coef = 1, n = Inf)
    diff.out <- na.omit(diff.out)
    
    diff.filter <- diff.out
    colnames(diff.filter)[c(1,4,5)] <- c('de.logFC','de.pval','de.qval')
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
      
      de.filter <- switch (cutoff.type,
                           param_0 = diff.filter,
                           param_1 = subset(diff.filter, diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,2]),
                           param_2 = subset(diff.filter, diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,2] & 
                                              diff.filter[,match(cutoff.s[2],colnames(diff.filter))] < cutoff[,3]))
    } else {
      cutoff.type <- paste('param', length(cutoff.s), sep = '_')
      
      de.filter <- switch (cutoff.type,
                           param_0 = diff.filter,
                           param_1 = subset(diff.filter, diff.filter[,match(cutoff.s,colnames(diff.filter))] < cutoff[,1]),
                           param_2 = subset(diff.filter, diff.filter[,match(cutoff.s[1],colnames(diff.filter))] < cutoff[,1] & 
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
    cat('\n',paste('The number of DEGs:',nrow(de.genes)),'\n')
  }
  
  diffcor.res <- DiffCor(exp.1, exp.2, tf2tar, cor.method = cor.method, p.adj = p.adj)

  diffcor.plus <- merge(diffcor.res,de.genes,by.x = 'Target',by.y = 'GeneSymbol')
  
  diffcor.plus <- diffcor.plus[,c(2,1,3:9)]
  diffcor.plus <- rbind(subset(diffcor.plus, diffcor.plus$high.condition == 1 & diffcor.plus$z.value < 0),
                        subset(diffcor.plus, diffcor.plus$high.condition == 2 & diffcor.plus$z.value > 0))

  return(diffcor.plus)
}



####################################################################################################
#' @name plotDysregExp
#' @title Visualize expression pattern of a dysreglaiton between conditions
#' @description Visualize expression pattern of a dysreglaiton between conditions.
#' @usage  plotDysregExp(tf, tar, exp.1, exp.2, exp1.lab, exp2.lab,
#'                       method, dysreg, conf.int.level, ...)
#' 
#' @param tf The TF of a gene dysregulation.
#' @param tar The target of a dysregulation.
#' @param exp.1 Expression matrix of a special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp.2 Expression matrix of an another special condition. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param exp1.lab The label condition of exp.1.
#' @param exp2.lab The label condition of exp.2.
#' @param method The mehod of dysregulation analysis, "dysreg", "diffcorplus", or "diffcor".
#' @param dysreg The results of dysreg.
#' @param conf.int.level Level controlling confidence region. 
#' @param ... further arguments to \code{ggscatter}.
#' 
#' @import ggpubr
#' 
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[test.genes,colnames(exp) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(exp),1000)
#' 
#' temp.1 <- exp.1[test.genes,]
#' temp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = temp.1, exp.2 = temp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', pValue = 0.01, ci = 0.90)
#' 
#' dysreg.res <- dysreg.out$dysreg
#' 
#' plotDysregExp(tf = dysreg.res$TF[1], tar = dysreg.res$Target[1],
#'               exp.1 = exp.1,exp.2 = exp.2, exp1.lab = 'PD',exp2.lab = 'RP',
#'               dysreg = dysreg.res, method ='dysreg', conf.int.level = 0.95)
#'               
#' }
#' 
#' @export
plotDysregExp <- function(tf, tar, exp.1, exp.2, exp1.lab, exp2.lab, 
                          method, dysreg, conf.int.level, ...){
  data=rbind(data.frame(TF=exp.1[rownames(exp.1)==tf,],
                        Target=exp.1[rownames(exp.1) ==tar,],
                        Group=exp1.lab,stringsAsFactors = F),
             data.frame(TF=exp.2[rownames(exp.2)==tf,],
                        Target=exp.2[rownames(exp.2)==tar,],
                        Group=exp2.lab,stringsAsFactors = F))

  if(method == 'dysreg'){
    reg.i=dysreg[dysreg$TF==tf & dysreg$Target==tar,]
    
    sp <- ggscatter(data, x = "TF", y = "Target",
                    add = "reg.line", conf.int = TRUE, conf.int.level = conf.int.level,
                    color = "Group", palette = rainbow(5)[c(4,1)],
                    shape = "Group") +  
      xlab(tf)+ylab(tar)
  } else {
    sp <- ggscatter(data, x = "TF", y = "Target",
                    add = "reg.line", conf.int = TRUE,
                    color = "Group", palette = rainbow(5)[c(4,1)],
                    shape = "Group") +                 
      xlab(tf)+ylab(tar)
  }
  
  sp <- switch (method,
    diffcor = sp + 
      stat_cor(aes(color = Group), method = cor.method, size = 6,
               label.x = min(data$TF) + (max(data$TF) - min(data$TF))*0.05,
               label.y = c(max(data$Target) + (max(data$Target) - min(data$Target))*0.2,
                           max(data$Target) + (max(data$Target) - min(data$Target))*0.1)),
    
    diffcorplus = sp + 
      stat_cor(aes(color = Group), method = cor.method, size = 6,
               label.x = min(data$TF) + (max(data$TF) - min(data$TF))*0.05,
               label.y = c(max(data$Target) + (max(data$Target) - min(data$Target))*0.2,
                           max(data$Target) + (max(data$Target) - min(data$Target))*0.1)),
  
    dysreg = sp + annotate("text", 
                           label = paste("Regulatory intensity: ",round(reg.i$unb.coef.y,3),'; ',
                                   paste(conf.int.level*100,'% CI: ', sep = ''),
                                   paste(round(reg.i$low.lim.y,3),round(reg.i$up.lim.y,3),sep = '~'),sep = ''),
                           x = min(data$TF) + (max(data$TF) - min(data$TF))*0.4, 
                           y = max(data$Target) + (max(data$Target) - min(data$Target))*0.2, 
                           size = 6, colour = rainbow(5)[1]) +
                  annotate("text", 
                           label = paste("Regulatory intensity: ",round(reg.i$unb.coef.x,3),'; ',
                                   paste(conf.int.level*100,'% CI: ', sep = ''),
                                   paste(round(reg.i$low.lim.x,3),round(reg.i$up.lim.x,3),sep = '~'),sep = ''),
                           x = min(data$TF) + (max(data$TF) - min(data$TF))*0.4, 
                           y = max(data$Target) + (max(data$Target) - min(data$Target))*0.1, 
                           size = 6, colour = rainbow(5)[4]))
  
  sp <- sp + theme_bw() +
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background  = element_rect(fill = "transparent")
    )   #backgroud
  
  sp <- sp+theme(
    axis.title.x = element_text(face = "bold",size = 16),
    axis.text.x  = element_text(size = 16),
    axis.title.y = element_text(face = "bold",size = 16),
    axis.text.y  = element_text(size = 16),
    legend.title=element_text(size=16,face = 'bold'),
    legend.text=element_text(size=16)
  )
  
  plot(sp)
}



####################################################################################################
#' @name RankDysReg
#' @title Rank dysregulatin
#' @description Rank dysregualtion by differential regulation and differential expression.
#' @usage RankDysReg(dysreg.res)
#' 
#' @param dysreg.res The results of \code{dysreg}.
#' 
#' 
#' @return The rank of dysregulations.
#' 
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[,colnames(exp) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(exp),1000)
#' 
#' temp.1 <- exp.1[test.genes,]
#' temp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = temp.1, exp.2 = temp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', pValue = 0.01, ci = 0.90)
#'                      
#' dysreg.res <- dysreg.out$dysreg
#' reg.rank <- RankDysReg(dysreg.res)
#' }
#' 
#' @export

RankDysReg <- function(dysreg.res){
  
  net.i <- dysreg.res
  net.i$weight <- abs((net.i$unb.coef.y - net.i$unb.coef.x) * net.i$logFC)
  dysreg.rank <- net.i[,c('TF','Target','weight')]
  
  dysreg.rank <- dysreg.rank[order(dysreg.rank$weight,decreasing = T),]
  
  rownames(dysreg.rank) <- NULL
  
  return(dysreg.rank)
}



####################################################################################################
#' @name RankDysTF
#' @title Rank TF with dysregulatin degree.
#' @description Rank TF with dysregualtion degree, which is measured by differential regulation and differential expression.
#' @usage RankDysTF(dysreg.res)
#' 
#' @param dysreg.res The results of \code{DysReg}.
#' 
#' @import igraph
#' 
#' @return The rank of TFs.
#' 
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[,colnames(exp) %in% group.2]
#' 
#' set.seed(1234)
#' test.genes <- sample(1:nrow(exp),1000)
#' 
#' temp.1 <- exp.1[test.genes,]
#' temp.2 <- exp.2[test.genes,]
#' 
#' dysreg.out <- DysReg(exp.1 = temp.1, exp.2 = temp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', pValue = 0.01, ci = 0.90)
#' 
#' dysreg.res <- dysreg.out$dysreg
#' 
#' tf.rank <- RankDysTF(dysreg.res)
#' }
#' 
#' @export

RankDysTF <- function(dysreg.res){
  
  net.i <- dysreg.res
  net.i$weight <- abs((net.i$unb.coef.y - net.i$unb.coef.x) * net.i$logFC)
  net.i <- graph.data.frame(net.i[,c('TF','Target','weight')], directed = TRUE)

  degree <- strength(net.i) #weighted degree
  degree <- as.data.frame(degree)
  degree$Gene <- rownames(degree)
  tf.dys <- degree[degree$Gene %in% dysreg.res$TF,]
  
  tf.rank <- tf.dys[order(tf.dys$degree,decreasing = T),]
  
  tf.rank <- tf.rank[,c(2,1)]
  rownames(tf.rank) <- NULL
  
  return(tf.rank)
}




####################################################################################################
#' @name combineDysReg
#' @title Search the best combination of dysregulations for building predictive signature.
#' @description Combine dysregulations for building predictive signature with genetic algorithm.
#' @usage combineDysReg(dysreg, exp.data, obj.func, pheno.data,
#'                      pop.size = 1000, select.rate=0.2, mut.rate = 0.1, add.rate=0.1, 
#'                      topN, train.rate = 0.6, iter=100)
#' 
#' @param dysreg The dysregulations output from \code{DysReg}.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param obj.func The objective function for genetic algorithm.
#' @param pheno.data The phenotype data corresponding to obj.func.
#' @param pop.size The size of populaiton generated in genetic algorithm.
#' @param select.rate The rate of selcting individual with good fitness from population.
#' @param mut.rate The rate of mutation.
#' @param add.rate The rate for adding individuals with not good fitness from population.
#' @param topN The top N individuals kept in each iteration.
#' @param train.rate The rate of sample for training model in cross-validation.
#' @param iter The times of iteration.
#' 
#' 
#' @import survival
#' @import survcomp
#' 
#' @return 
#' The results of gene dysregulation analysis:
#'  \item{individ}{The combination of dysregulaitons.}
#'  \item{best.fitness}{The fitness value for each individual.}
#'  
#' 
#' @examples
#' \donttest{
#' data(exp)
#' exp[1:5,1:5]
#' 
#' data(tf2tar)
#' head(tf2tar)
#' 
#' group.1 <- clin_data$sample[which(clin_data$ResponseLevel == 'PD')]
#' exp.1 <- exp[,colnames(exp) %in% group.1]
#' 
#' group.2 <- clin_data$sample[which(clin_data$ResponseLevel %in% c('CR','PR','SD'))]
#' exp.2 <- exp[,colnames(exp) %in% group.2]
#' 
#' dysreg.out <- DysReg(exp.1 = exp.1, exp.2 = exp.2, tf2tar, 
#'                      de.genes = NULL, de.pval = 0.05, 
#'                      grn.method = 'Boruta', pValue = 0.01, ci = 0.90)
#'
#' dysreg.res <- dysreg.out$dysreg
#' dysreg <- dysreg.res[,1:2]
#' 
#' data(clin_data)
#' head(clin_data)
#' pheno.data <- clin_data[,c(1,3,4)]
#' pheno.data <- pheno.data[!is.na(pheno.data$os) & !is.na(pheno.data$censOS),]
#' 
#' combdysreg.out <- combineDysReg(dysreg = dysreg, exp.data = exp,
#'                                 obj.func = 'Obj.Cindex', pheno.data = pheno.data, 
#'                                 pop.size = 1000, select.rate=0.2, 
#'                                 mut.rate=0.1, add.rate=0.1, 
#'                                 topN = 10, train.rate = 0.6, iter = 100)
#' 
#' best.fit <- combdysreg.out$best.fitness
#' best.individ <- combdysreg.out$individ
#' head(best.fit)
#' summary(best.fit$fitness)
#' 
#' ##Check the individ with the best fitness
#' best.individ[which.max(best.fit$fitness),]
#'                                                                  
#' }
#' 
#' @export
combineDysReg <- function(dysreg,exp.data,obj.func, pheno.data, pop.size = 1000,
                          select.rate=0.2,mut.rate=0.1,add.rate=0.1, topN = 10,
                          train.rate = 0.6,iter = 100){
  
  chr.len <- nrow(dysreg)
  if(chr.len > nrow(pheno.data)*0.5){
    stop('The number of inputting dysregulations is larger than half of example size, please filter dysregulations.')
  }
  
  ## generate a random initial population
  pop <- do.call(rbind,lapply(1:pop.size,function(i,chr.len){
    gene.num <- sample(1:chr.len,1,replace = F)
    individ <- sample(c(rep(1,gene.num),rep(0, chr.len - gene.num)),chr.len,replace = F)
    return(individ)
  },chr.len = nrow(dysreg)))
  
  
  best.output <- list()
  best.output$individ <- vector()
  best.output$best.fitness <- vector()
  
  for (i in 1: iter) {
    
    rownames(pop) <- c(1:nrow(pop))
    
    ##calculate fitness
    obj.fitness <- apply(pop, 1, obj.func,
                         dysreg = dysreg,exp.data = exp.data,pheno.data = pheno.data,train.rate = train.rate)
    obj.fitness <- as.data.frame(obj.fitness)
    obj.fitness$indival <- rownames(obj.fitness)
    obj.fitness <- obj.fitness[,c(2,1)]
    colnames(obj.fitness)[2] <- 'obj.fit'
    
    pop.obj <- cbind(pop,obj.fitness)
    pop.obj <- pop.obj[order(pop.obj$obj.fit,decreasing = T),]
    obj.value <- obj.fitness[order(obj.fitness$obj.fit,decreasing = T),]
    
    best.fitness <- data.frame(iter = i,fitness = obj.value$obj.fit[1:topN],stringsAsFactors = F)
    best.output$best.fitness <- rbind(best.output$best.fitness,best.fitness)
    
    best.individ <- pop.obj[pop.obj$indival %in% obj.value$indival[1:topN],1:(ncol(pop.obj)-2)]
    best.output$individ <- rbind(best.output$individ,best.individ)
    
    
    ## evolution prorecess
    if(i < (iter-1)){
      children.pop <- evolution(pop = pop, pop.obj = pop.obj, obj.value = obj.value, 
                                pop.size = pop.size, select.rate = 0.2, mut.rate = 0.1,add.rate = 0.1)
      
      pop <- children.pop
    }
  }
  
  return(best.output)
}

evolution <- function(pop, pop.obj, obj.value, pop.size = pop.size, 
                      select.rate = select.rate, mut.rate = mut.rate, add.rate = add.rate){
  
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
    
    individ.1 <- as.numeric(selected.pop[selected.pop$indival==individ.1,1:(ncol(selected.pop)-2)])
    individ.2 <- as.numeric(selected.pop[selected.pop$indival==individ.2,1:(ncol(selected.pop)-2)])
    
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



####################################################################################################
#' @name Obj.Cindex
#' @title Calculate C-Index in survival data for using as Objective function.
#' @description Calculate C-Index with cross-validation in survival data for using as Objective function.
#' @usage Obj.Cindex(combine.genes, dysreg, exp.data, 
#'                   pheno.data, train.rate)
#' 
#' @param combine.genes The combination of genes.
#' @param dysreg The dysregulations output from \code{DysReg}.
#' @param exp.data Expression matrix. Columns correspond to genes, rows correspond to experiments. The matrix is expected to be already normalized.
#' @param pheno.data The phenotype data corresponding to objective function. At here, it is survival data.
#' @param train.rate The rate of sample for training model in cross-validation.
#' 
#' @import survival
#' @import survcomp
#' 
#' @return 
#' The fitness, C-Index, for each individual.
#' 
#' @export
Obj.Cindex <- function(combine.genes,dysreg,exp.data,pheno.data,train.rate){
  regs <- dysreg[combine.genes !=0,]
  genes <- union(regs[,1],regs[,2])
  
  exp <- exp.data[rownames(exp.data) %in% genes,]
  exp <- as.data.frame(t(exp),stringsAsFactors = F)
  exp$sample <- rownames(exp)
  colnames(exp) <- sub('-','_',colnames(exp))
  
  colnames(pheno.data)[2:3] <- c('time','event')
  pheno.data <- pheno.data[!is.na(pheno.data$time) & !is.na(pheno.data$event),] 
  
  data <- merge(exp, pheno.data, by = 'sample')
  data <- data[,-1]
  
  cv.ci <- do.call(rbind,lapply(1:10,function(i,data, train.rate){
    cv.select <- sample(1:nrow(data),round(nrow(data)*train.rate),replace = F)
    train.data <- data[cv.select,]
    test.data <- data[setdiff(1:nrow(data),cv.select),]
    
    genes <- setdiff(colnames(exp),'sample')
    ff <- as.formula(paste("Surv(time,event)~", paste(genes, collapse = "+")))
    train.mod <- coxph(ff, data=train.data)
    train.mod$coefficients[is.na(train.mod$coefficients)] <- 0 
    score <- as.matrix(test.data[,1:(ncol(test.data) - 2)])%*% train.mod$coefficients
    res.i <- data.frame(times=i,fitness=concordance.index(score,test.data$time,test.data$event)$c.index)
    return(res.i)
  },data = data, train.rate = train.rate))
  
  cv.ci <- mean(cv.ci$fitness,trim = 0.1)
  return(cv.ci)
}



####################################################################################################
#' @name exp
#' @docType data
#' @title Example of expression data.
#' @description An example of expression data.
#' @usage data(exp.1)
#' 
#' @details The example data was derived from IMvigor210CoreBiologies, which includes RNA-seq data of 348 samples and corresponding drug response of atezolizumab.
#' @examples
#' data(exp)
#' exp[1:5,1:5,]
#' 
#' @keywords datasets
#' @export



####################################################################################################
#' @name clin_data
#' @docType data
#' @title Example of clinical data.
#' @description An example of clinical data.
#' @usage data(clin_data)
#' 
#' @details The example data was derived from IMvigor210CoreBiologies, which includes RNA-seq data of 348 samples and corresponding drug response of atezolizumab.
#' @examples
#' data(clin_data)
#' head(clin_data)
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name tf2tar
#' @docType data
#' @title TF-target regulations.
#' @description TF-target regulations within an reference GRN.
#' @usage data(exp.1)
#' 
#' @details The TF-target relationships were derived from predicting by FIMO.
#' @examples
#' data(tf2tar)
#' head(tf2tar)
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name de_genes
#' @docType data
#' @title Example of DEGs input n function DysReg and DiffCorPlus.
#' @description Example of DEGs input n function DysReg and DiffCorPlus.
#' @usage data(de_genes)
#' 
#' @details de_genes contains three columns, GeneSymbol is the gene symbol of DEGs, high.condition indicates which comdition expresses high level, logFC is logFC value.
#' @examples
#' data(de_genes)
#' head(de_genes)
#' 
#' @keywords datasets
#' @export

