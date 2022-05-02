#************************
#
#
# v.0.1.2  
# daniel.frank@ucdenver.edu
# 03 Mar 2019
#
#
#************************

library("car")

otu.test_ra <- function(obj, predictor, mask="", main="", plot="p", 
                        drop.zeros=TRUE, filen="", pdfx = 8, pdfy = 2.5, lwd = 2, rm.other = FALSE) {
  if (length(mask) == 1) mask <- !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  
  if (plot == "none") filen = "" #turn off PDF saving

  pred_mask = predictor[mask==TRUE]
  obj_mask = otu.subset_libs(obj, mask, drop.zeros = drop.zeros)

  if (is.factor(predictor) == TRUE) pred_mask = droplevels(pred_mask)
  
  signs = ""

  opar = par()
  if (plot != "none") {
    if (plot == "both") {
      par(mfrow=c(2,1), mar=c(3, 4, 3, 4))
      pdfy = 2 * pdfy
    }
    else par(mfrow=c(1,1), mar=c(3, 4, 3, 4))
  }
  
  if (is.factor(predictor) == TRUE & (length(levels(pred_mask)) == 2)) {
    cat("Found factor with two levels: Perform Wilcoxon test\n")
    res <- otu.mul_wilcox(obj_mask$per, pred_mask, rm.other = rm.other)
    ret = cbind(res$OTU, res$p_value, res$fdr, res$med1, res$med2, res$stat)
    colnames(ret) <- c("OTU", "Pvalue", "FDR", paste("Median:Grp1=",res$grp1,sep=""), paste("Median:Grp2=",res$grp2, sep=""), "W_Statistic")
    
    signs = (res$med1 > res$med2)
    if (plot != "none") {
      if (plot == "both" || plot == "p") otu.plot_pvalue_bar(as.numeric(res$p_value), signs = signs, main=paste(main, "(P-value)"), ylab="-log10(P)", lwd=lwd)
      if (plot == "both" || plot == "fdr") otu.plot_pvalue_bar(as.numeric(res$fdr), signs = signs, main=paste(main, "(FDR)"), ylab="-log10(FDR)", lwd=lwd)
    }
  } else if (is.factor(pred_mask) == TRUE) {
    cat("Found factor with >2 levels: Perform Kruskal-Wallis test\n")
    res <- otu.mul_krusk(obj_mask$per, pred_mask, rm.other = rm.other)
    ret = cbind(res$OTU, res$p_value, res$fdr, res$stat)
    colnames(ret) <- c("OTU", "Pvalue", "FDR", "ChiSq_Statistic")
    if (plot != "none") {
      if (plot == "both" || plot == "p") otu.plot_pvalue_bar(as.numeric(res$p_value), signs = signs, main=paste(main, "(P-value)"), ylab="-log10(P)", lwd=lwd)
      if (plot == "both" || plot == "fdr") otu.plot_pvalue_bar(as.numeric(res$fdr), signs = signs, main=paste(main, "(FDR)"), ylab="-log10(FDR)", lwd=lwd)
    }
  }else if (is.numeric(pred_mask) == TRUE) {
    cat("Found numerical predictor\n")
    res <- otu.mul_corr(obj_mask$per, pred_mask, method="spearman", rm.other = rm.other)
    ret = cbind(res$OTU, res$p_value, res$fdr, signif(res$rho,2), res$statistic)
    #colnames(ret) <- c("OTU", "Pvalue", "FDR", "Rho","Statistic")
    colnames(ret) <- c("OTU", "Pvalue", "FDR", "Rho")
    
    if (plot != "none") {
      signs = (res$rho >= 0)
      
      if (plot == "p" || plot == "both") {
        otu.plot_pvalue_bar(as.numeric(res$p_value), signs = signs, main=paste(main, "(P-value)"), ylab="-log10(P)", lwd=lwd)
      } else if (plot == "fdr" || plot == "both") {
        otu.plot_pvalue_bar(as.numeric(res$fdr), signs = signs, main=paste(main, "(FDR)"), ylab="-log10(FDR)", lwd=lwd)
      }
      
      signs = ""
      if (plot == "rho" || plot == "both") otu.plot_pvalue_bar(as.numeric(res$rho), signs = signs, main=paste(main, "(Rho)"), ylab="rho", logtrans=FALSE, plot="cor", lwd=lwd)
    }
  }
  
  if (filen != "") {
    dev.copy(pdf, filen, width = pdfx, height = pdfy)
    dev.off()
  }
  
  par(opar)
  return(ret)
}

otu.test_ra_lm <- function(obj, predictor, factor = "", mask="", main="", 
                           plot="none", drop.zeros=TRUE, filen="", pdfx = 8, pdfy = 2.5, lwd = 2, rm.other = FALSE) {
  if (length(mask) == 1) mask <- !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  
  if (plot == "none") filen = "" #turn off PDF saving
  
  pred_mask = predictor[mask==TRUE & !is.na(mask)]
  
  obj_mask = otu.subset_libs(obj, mask, drop.zeros = drop.zeros)
  if (is.factor(predictor) == TRUE) pred_mask = droplevels(pred_mask)
  
  #otu.summary(obj_mask)
  
  fac_mask = ""
  if (length(factor) > 1) {
    fac_mask = factor[mask==TRUE & !is.na(mask)]
    if (is.factor(factor) == TRUE) fac_mask = droplevels(fac_mask)
  }
  #print(fac_mask)
  signs = ""
  
  opar = par()
  if (plot != "none") {
    if (plot == "both") {
      par(mfrow=c(2,1), mar=c(3, 4, 3, 4))
      pdfy = 2 * pdfy
    }
    else par(mfrow=c(1,1), mar=c(3, 4, 3, 4))
  }
  
  if (is.factor(predictor) == TRUE & (length(levels(pred_mask)) == 2)) {
    cat("Found factor with two levels\n")
    #Calculate median values (otu.mul_lm uses clr transformed data, so medians aren't relative abundances)
    
    n <- ncol(obj_mask$per)

    med1 <- c(1:n)
    med2 <- c(1:n)
    levels = levels(pred_mask)
    
    for (i in 1:n) {
      vec = (obj_mask$per)[,i]
      x = subset(vec, pred_mask == levels[1])
      y = subset(vec, pred_mask == levels[2])
      med1[i] = round(100*median(x, na.rm = TRUE), 4)
      med2[i] = round(100*median(y, na.rm = TRUE), 4)
    }
    
    #res <- otu.mul_lm0(obj_mask$clr, pred_mask, fac_mask)
    res <- otu.mul_lm(obj_mask$clr, pred_mask, fac_mask)
    
    ret = cbind(res$OTU, res$p_value1, res$fdr1, med1, med2)
    colnames(ret) <- c("OTU", "Pvalue", "FDR", paste("Median:Grp1=",res$grp1,sep=""), paste("Median:Grp2=",res$grp2, sep=""))
    
    signs = (med1 > med2)
    if (plot != "none") {
      if (plot == "both" || plot == "p") otu.plot_pvalue_bar(as.numeric(res$p_value1), signs = signs, main=paste(main, "(P-value)"), ylab="-log10(P)", lwd=lwd)
      if (plot == "both" || plot == "fdr") otu.plot_pvalue_bar(as.numeric(res$fdr1), signs = signs, main=paste(main, "(FDR)"), ylab="-log10(FDR)", lwd=lwd)
    }
  } else if (is.factor(pred_mask) == TRUE) {
    cat("Found factor with >2 levels")
    res <- otu.mul_lm(obj_mask$per, pred_mask, fac_mask)
    ret = cbind(res$OTU, res$p_value1, res$fdr1)
    colnames(ret) <- c("OTU", "Pvalue", "FDR")
    if (plot != "none") {
      if (plot == "both" || plot == "p") otu.plot_pvalue_bar(as.numeric(res$p_value1), signs = signs, main=paste(main, "(P-value)"), ylab="-log10(P)", lwd=lwd)
      if (plot == "both" || plot == "fdr") otu.plot_pvalue_bar(as.numeric(res$fdr1), signs = signs, main=paste(main, "(FDR)"), ylab="-log10(FDR)", lwd=lwd)
    }
  }
  #}else if (is.numeric(pred_mask) == TRUE) {   #Implement later
  #   cat("Found numerical predictor\n")
  #  res <- otu.mul_corr(obj_mask$per, pred_mask, method="spearman")
  #  ret = cbind(res$OTU, res$p_value, res$fdr, signif(res$rho,2))
  #  colnames(ret) <- c("OTU", "Pvalue", "FDR", "Rho")
  #  
  #  if (plot != "none") {
  #    signs = (res$rho >= 0)
  #    
  #    if (plot == "p" || plot == "both") otu.plot_pvalue_bar(as.numeric(res$p_value), signs = signs, main=paste(main, "(P-value)"), ylab="-log10(P)")
  #    else if (plot == "fdr" || plot == "both") otu.plot_pvalue_bar(as.numeric(res$fdr), signs = signs, main=paste(main, "(FDR)"), ylab="-log10(FDR)")
  #    
  #    signs = ""
  #    if (plot == "rho" || plot == "both") otu.plot_pvalue_bar(as.numeric(res$rho), signs = signs, main=paste(main, "(Rho)"), ylab="rho", logtrans=FALSE, plot="cor")
  #  }
  #}
  
  if (filen != "") {
    dev.copy(pdf, filen, width = pdfx, height = pdfy)
    dev.off()
  }
  
  par(opar)
  return(ret)
}

otu.summarize_princo <- function(obj) {
  cat("\nObject: ",deparse(substitute(obj)),"\n\n")
  cat("***CORRELATION SUMMARY:\n")
  print((summary(obj$cor))$importance)

  cat("\n\n***COVARIATION SUMMARY:\n")
  print((summary(obj$cov))$importance)
  cat("\n\n***BRAY-CURTIS SUMMARY:\n")
  print(summary(eigenvals(obj$bc)))
  cat("\n\n***MORISITA_HORN SUMMARY:\n")
  print(summary(eigenvals(obj$mh)))
  cat("\n\n")
}

otu.mul_princo <- function(obj) {
  cor = prcomp(obj$clr, scale=TRUE, center=TRUE)
  cov = prcomp(obj$clr) # pc1-3: 34% variance
  bcdst = vegdist(obj$cts, method="bray")
  bc = wcmdscale(bcdst, eig = TRUE, x.ret = TRUE) #vegan call of cmdscale
  mhdst = vegdist(obj$cts, method="horn")
  mh = wcmdscale(mhdst, eig = TRUE, x.ret = TRUE) #vegan call of cmdscale
  
  return(list(cor = cor, cov = cov, bc = bc, mh = mh))
}

otu.mul_plot_princo <- function(obj = obj, 
                            pc1 = 1, 
                            pc2 = 2, 
                            aes = vector, 
                            mask = "", 
                            main="", 
                            lgd.title="Legend", 
                            drop.levels = FALSE, 
                            cex.mul = 1, 
                            plot.met = 0,
                            num.cols = -1,
                            filen="") {

  if (length(mask) == 1) mask <- !vector(mode = "logical", length = length(aes)) #new vector is filled with "FALSE", so take inverse
  aes <- aes[mask]
  
  knd = ""
  cpa = ""
  if (class(aes) == "factor") {
    print("Factor")
    knd = "factor"
    if (drop.levels == TRUE)  aes <- droplevels(aes) 
    cpa = col_pch_assign_factor(aes, num.cols)
  } else if (class(aes) == "numeric") {
    print("Numeric")
    knd = "numeric"
    cpa = col_pch_assign_numeric(aes, num.cols)
  }
  corx = (obj$cor)$x[mask,pc1]
  cory = (obj$cor)$x[mask,pc2]
  covx = (obj$cov)$x[mask,pc1]
  covy = (obj$cov)$x[mask,pc2]
  mhx = (obj$mh)$points[mask,pc1]
  mhy = (obj$mh)$points[mask,pc2]
  bcx = (obj$bc)$points[mask,pc1]
  bcy = (obj$bc)$points[mask,pc2]
  
  xl = paste("PC", as.character(pc1), sep="")
  yl = paste("PC", as.character(pc2), sep="")
  
  if (plot.met == 0) {
    col = "black"
    bg =  cpa$colvec
    lwd = 1.0
  } else if (plot.met == 1) {
    col = cpa$colvec
    bg = "white"
    lwd = 2.0
  }
  pch = cpa$pchvec
  
  opar = par()
  par(mfrow=c(2,3), mar=c(4,4, 2, 0.5))
  plot(corx, cory, xlab = xl, ylab=yl, type="n", main ="PCA-COR")
  points(corx, cory, col = col,  bg=bg, pch = pch, cex=cex.mul * cpa$cexvec, lwd = lwd)

  plot(covx, covy, xlab = xl, ylab=yl, type="n", main ="PCA-COV")
  points(covx, covy, col = col,  bg=bg, pch = pch, cex=cex.mul * cpa$cexvec, lwd = lwd)
  
  #legend("left",legend= cpa$lgdlvl, title=lgd.title, bg = "grey96", pt.bg=cpa$lgdcol, pch = cpa$lgdpch, cex=cx)
  if (knd == "factor") {
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    cx=1
    if (length(cpa$lgdlvl) > 10) cx = 0.6
    legend("left",legend= cpa$lgdlvl, title=lgd.title, bg = "white", pt.bg=cpa$lgdcol, pch = cpa$lgdpch, cex=cx)
  } else if (knd == "numeric") {
    colfunc = colorRampPalette(c("blue", "white", "red"))
    color.bar(colfunc(100), min = 0, max = 1, nticks=5, lgd.title=lgd.title)
  }
  plot(bcx, bcy, xlab = xl, ylab=yl, type="n", main ="PCoA-BC")
  points(bcx, bcy, col = col,  bg=bg, pch = pch, cex=cex.mul * cpa$cexvec, lwd = lwd)
  
  plot(mhx, mhy, xlab = xl, ylab=yl, type="n", main ="PCoA-MH")
  points(mhx, mhy, col = col,  bg=bg, pch = pch, cex=cex.mul * cpa$cexvec, lwd = lwd)
  

  if (filen != "") {
    #dev.copy(pdf, filen, width = pdfx, height = pdfy)
    dev.copy(pdf, filen)
    dev.off()
  }
  
  par(opar)
}

otu.mul_krusk <- function(frame = frame, factor = factor, rm.other = FALSE)
{
  if (nrow(frame) != length(factor)) stop("Frame and factor must have same lengths")
  if (rm.other == TRUE) {
    frame = frame[, colnames(frame) != "Other"]
  }
  n <- ncol(frame)
  
  if (is.factor(factor) == 0) factor <- as.factor(factor)
  tmp = factor(factor)
  factor = tmp
  levels <- levels(factor) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 
  
  if (level_count < 2) stop("Factor must have > 1 levels")
  
  col_name <- c(1:n)
  krusk_p_value <- c(1:n)
  statistic <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value <- c(1:n)
  
  for (i in 1:n) {
    vec = frame[,i]
    ks <- kruskal.test(frame[,i] ~ factor)
    krusk_p_value[i] = signif(ks$p.value, 2)
    statistic[i] = signif(ks$statistic,2)
  }
  
  fdr_p_value <- p.adjust(krusk_p_value, method = "fdr")
  fdr_p_value <- signif(fdr_p_value, 2)
  return(list(OTU = colnames(frame), p_value = krusk_p_value, fdr = fdr_p_value, stat = statistic))
}

otu.mul_wilcox <- function(frame = frame, factor = factor, rm.other = FALSE)
{
  if (nrow(frame) != length(factor)) stop("Frame and factor must have same lengths")
  if (rm.other == TRUE) {
    frame = frame[, colnames(frame) != "Other"]
  }
  n <- ncol(frame)
  
  if (is.factor(factor) == 0) factor <- as.factor(factor)
  tmp = factor(factor)
  factor = tmp
  levels <- levels(factor) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 
  
  if (level_count < 2) stop("Factor must have > 1 levels")
  
  col_name <- c(1:n)
  wilcox_p_value <- c(1:n)
  statistic <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value <- c(1:n)
  med1 <- c(1:n)
  med2 <- c(1:n)
  
  for (i in 1:n) {
    vec = frame[,i]
    ks <- wilcox.test(frame[,i] ~ factor)
    wilcox_p_value[i] = signif(ks$p.value, 2)
    statistic[i] = signif(ks$statistic,3)
    
    x = subset(vec, factor == levels[1])
    y = subset(vec, factor == levels[2])
    med1[i] = round(100*median(x, na.rm = TRUE), 4)
    med2[i] = round(100*median(y, na.rm = TRUE), 4)
  }
  
  fdr_p_value <- p.adjust(wilcox_p_value, method = "fdr")
  fdr_p_value <- signif(fdr_p_value, 2)
  return(list(OTU = colnames(frame), p_value = wilcox_p_value, fdr = fdr_p_value, stat = statistic, med1 = med1, med2 = med2, grp1 = levels[1], grp2 = levels[2]))
}

otu.mul_lm0 <- function(frame = frame, factor1 = factor, factor2 = "")
{
  if (nrow(frame) != length(factor1)) stop("Frame and factor must have same lengths")
  n <- ncol(frame)
  
  if (is.factor(factor1) == FALSE) factor1 <- as.factor(factor1)
 
  levels <- levels(factor1) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 

  if (level_count < 2) stop("Factor must have > 1 levels")
  
  col_name <- c(1:n)
  lm_p_value1 <- c(1:n)
  lm_p_value2 <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value1 <- c(1:n)
  fdr_p_value2 <- c(1:n)
  med1 <- c(1:n)
  med2 <- c(1:n)
  statistic <- c(1:n)
  
  for (i in 1:n)
  {
    if (length(factor2) <= 1) {
      result = coefficients(summary(lm(frame[,i] ~ factor1, na.action = na.omit)))
    } else {
      result = coefficients(summary(lm(frame[,i] ~ factor1 + factor2, na.action = na.omit)))
    }
    # print(result)
    lm_p_value1[i] = signif(result[2,4], 2)
    
    vec = frame[,i]
    x = subset(vec, factor1 == levels[1])
    y = subset(vec, factor1 == levels[2])
    med1[i] = round(100*median(x, na.rm = TRUE), 4)
    med2[i] = round(100*median(y, na.rm = TRUE), 4)
  }
  
  fdr_p_value1 <- p.adjust(lm_p_value1, method = "fdr")
  fdr_p_value1 <- signif(fdr_p_value1, 2)
  
  fdr_p_value2 <- p.adjust(lm_p_value2, method = "fdr")
  fdr_p_value2 <- signif(fdr_p_value2, 2)
  
  return(list(OTU = colnames(frame), p_value1 = lm_p_value1, fdr1 = fdr_p_value1,  p_value2 = lm_p_value2, fdr2 = fdr_p_value2, med1 = med1, med2 = med2, grp1 = levels[1], grp2 = levels[2]))
  #return(list(OTU = colnames(frame), p_value = wilcox_p_value, fdr = fdr_p_value, stat = statistic, med1 = med1, med2 = med2, grp1 = levels[1], grp2 = levels[2]))
}

otu.mul_lm <- function(frame = frame, factor1 = factor, factor2 = "")
{
  if (nrow(frame) != length(factor1)) stop("Frame and factor must have same lengths")
  n <- ncol(frame)
  
  if (is.factor(factor1) == FALSE) factor1 <- as.factor(factor1)
  
  levels <- levels(factor1) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 
  
  if (level_count < 2) stop("Factor must have > 1 levels")
  
  col_name <- c(1:n)
  lm_p_value1 <- c(1:n)
  lm_p_value2 <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value1 <- c(1:n)
  fdr_p_value2 <- c(1:n)
  med1 <- c(1:n)
  med2 <- c(1:n)
  statistic <- c(1:n)
  
  for (i in 1:n)
  {
    if (length(factor2) <= 1) {
      result = coefficients(summary(lm(frame[,i] ~ factor1, na.action = na.omit)))
      # print(result)
      lm_p_value1[i] = signif(result[2,4], 2)
    } else {
      mod = lm(frame[,i] ~ factor1 + factor2, na.action = na.omit)
      anova_mod = car::Anova(mod, type = 2)
      #Extract Type 2 pvalue for factor1
      lm_p_value1[i] = signif(anova_mod[1, "Pr(>F)"], 2)
      lm_p_value2[i] = signif(anova_mod[2, "Pr(>F)"], 2)
    }
    
    
    vec = frame[,i]
    x = subset(vec, factor1 == levels[1])
    y = subset(vec, factor1 == levels[2])
    med1[i] = round(100*median(x, na.rm = TRUE), 4)
    med2[i] = round(100*median(y, na.rm = TRUE), 4)
  }
  
  fdr_p_value1 <- p.adjust(lm_p_value1, method = "fdr")
  fdr_p_value1 <- signif(fdr_p_value1, 2)
  
  fdr_p_value2 <- p.adjust(lm_p_value2, method = "fdr")
  fdr_p_value2 <- signif(fdr_p_value2, 2)
  
  return(list(OTU = colnames(frame), p_value1 = lm_p_value1, fdr1 = fdr_p_value1,  p_value2 = lm_p_value2, fdr2 = fdr_p_value2, med1 = med1, med2 = med2, grp1 = levels[1], grp2 = levels[2]))
  #return(list(OTU = colnames(frame), p_value = wilcox_p_value, fdr = fdr_p_value, stat = statistic, med1 = med1, med2 = med2, grp1 = levels[1], grp2 = levels[2]))
}

otu.mul_adonis <- function(obj = obj, factor = factor, mask = "", method="bray", 
                           permutations = 1000, filen = "", plotp = FALSE, bin = FALSE) {
  if (!is.factor(factor)) stop("Must enter factor (try as.factor(factor))")
  if (length(mask) <= 1) mask <- !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  roundto = 5
  #build new mask that takes into account any "na" in the factor
  combined_mask = (mask & !is.na(mask)) & !is.na(factor)

  obj2= otu.subset_libs(obj, mask = combined_mask,drop.zeros=TRUE, vrb=TRUE)
  if (bin == TRUE) {
    cts = (obj2$bin)
  } else {
    cts = (obj2$cts)
  }
  fac = droplevels(factor[combined_mask])
  
  #Set up data structures
  levs = levels(fac)
  levs_len = length(levs)
  nstr = paste("Total:", levs_len, " levels")
  
  pmat = matrix(nrow = levs_len, ncol = levs_len, dimnames = list(levs,levs))
  #qmat = matrix(nrow = levs_len, ncol = levs_len, dimnames = list(levs,levs))
  smat = matrix(nrow = levs_len, ncol = levs_len, dimnames = list(levs,levs))
  #qvals = c(1:(levs_len * (levs_len - 1))/2)
  perms = (levs_len * (levs_len - 1))/2
  pvals = c(1:perms)
    
  if (filen != "") sink(filen,split="TRUE")
  
  cat("======================\nAll levels:", deparse(substitute(factor)),"\n")
  cat("Factors:",levels(fac),"\n ", nstr,"\n")
      
  res = adonis(cts ~ fac, method = method, permutations = permutations)
  print(res)
  cat("======================\n\n")
  
  if (levs_len > 2) {
    q=1
    for (i in 1:(levs_len-1)) {
      for (j in (i + 1):levs_len) {
        lev1 = levs[i]
        lev2 = levs[j]
        n1 = sum(fac == lev1 )
        n2 = sum(fac == lev2 )
        
        msk = (fac == lev1 | fac == lev2)
        f = fac[msk]
        o = cts[msk, ]
        
        #f = fac[fac == lev1 | fac == lev2]
        #o = cts[fac == lev1 | fac == lev2, ]
        
        #n1str = paste("level1: ", lev1, " (n=", n1, ") ", sep="")
        #n2str = paste("  level2: ", lev2, " (n=", n2, ") ", sep="")
        
        n1str = paste("Level ", i, ": ", lev1, " (n=", n1, ") ", sep="")
        n2str = paste("  Level ", j, ": ", lev2, " (n=", n2, ") ", sep="")
        #n2str = paste("  level2: ", lev2, " (n=", n2, ") ", sep="")
        
        cat("======================\n", n1str,n2str,  "   ", nstr, "\n")
      
        res = adonis(o ~ f, method = method, permutations = permutations)
        #print(res)
        #print(str(res))
        p = ((res$aov.tab)[["Pr(>F)"]])[1]
        if (is.na(p)) p = 1.0
        
        pmat[i,j] = round(p,roundto)
        pmat[j,i] = round(p,roundto)
        pvals[q] = round(p,roundto)
        q = q + 1
        
        if (p < 0.001) {
          smat[i,j] = "***"
          smat[j,i] = "***"
        } else if (p < 0.01) {
          smat[i,j] = "**"
          smat[j,i] = "**"
        } else if (p < 0.05) {
          smat[i,j] = "*"
          smat[j,i] = "*"
        } else if (p < 0.1) {
          smat[i,j] = "."
          smat[j,i] = "."
        } else {
          smat[i,j] = "-"
          smat[j,i] = "-"
        }
        
        print(res)
        cat("======================\n\n")
      }
    }
    for (x in 1:levs_len) {
      smat[x,x] = "-"
      pmat[x,x] = 1.0
    }
  
    qvals = p.adjust(pvals, method="fdr")
    qvals = round(qvals,roundto)
    q=1
    for (i in 1:(levs_len-1)) {
      for (j in (i + 1):levs_len) {
        pmat[i, j] = qvals[q]
        if (qvals[q] < 0.001) {
         smat[i,j] = "***"
        } else if (qvals[q] < 0.01) {
          smat[i,j] = "**"
        } else if (qvals[q] < 0.05) {
          smat[i,j] = "*"
        } else if (qvals[q] < 0.1) {
          smat[i,j] = "."
        } else {
          smat[i,j] = "-"
        }
        q = q + 1
      }
    }
    cat("")
    cat("P-value and FDR matrix:\n")
    print(pmat)
    cat("")
    cat("    P = lower left triangle\n")
    cat("    FDR = upper right triangle\n")
    cat("======================\n\n")
    cat("Symbol matrix:\n")
    print(noquote(smat))
    cat("======================\n\n")
  
  }
  
  if (filen != "")  sink()
  
  if (plotp == TRUE && levs_len > 2) otu.plot_pvalue_heat(pmat, maxp = 0.05)
  return(list(p=pmat,sym=smat))
}

otu.mul_adonis2 <- function(obj = obj, factor = factor, mask = "", method="bray", permutations = 1000, filen = "", plotp = FALSE) {
  if (!is.factor(factor)) stop("Must enter factor (try as.factor(factor))")
  if (length(mask) <= 1) mask <- !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  roundto = 5
  #build new mask that takes into account any "na" in the factor
  combined_mask = (mask & !is.na(mask)) & !is.na(factor)
  
  obj2= otu.subset_libs(obj, mask = combined_mask,drop.zeros=TRUE, vrb=TRUE)
  cts = (obj2$cts)
  fac = droplevels(factor[combined_mask])
  
  #Set up data structures
  levs = levels(fac)
  levs_len = length(levs)
  nstr = paste("Total:", levs_len, " levels")
  
  pmat = matrix(nrow = levs_len, ncol = levs_len, dimnames = list(levs,levs))
  #qmat = matrix(nrow = levs_len, ncol = levs_len, dimnames = list(levs,levs))
  smat = matrix(nrow = levs_len, ncol = levs_len, dimnames = list(levs,levs))
  #qvals = c(1:(levs_len * (levs_len - 1))/2)
  perms = (levs_len * (levs_len - 1))/2
  pvals = c(1:perms)
  
  if (filen != "") sink(filen,split="TRUE")
  
  cat("======================\nAll levels:", deparse(substitute(factor)),"\n")
  cat("Factors:",levels(fac),"\n ", nstr,"\n")
  
  res = adonis(cts ~ fac, method = method, permutations = permutations)
  print(res)
  cat("======================\n\n")
  
  if (levs_len > 2) {
    q=1
    for (i in 1:(levs_len-1)) {
      for (j in (i + 1):levs_len) {
        lev1 = levs[i]
        lev2 = levs[j]
        n1 = sum(fac == lev1 )
        n2 = sum(fac == lev2 )
        
        msk = (fac == lev1 | fac == lev2)
        f = fac[msk]
        o = cts[msk, ]
        
        n1str = paste("Level ", i, ": ", lev1, " (n=", n1, ") ", sep="")
        n2str = paste("  Level ", j, ": ", lev2, " (n=", n2, ") ", sep="")
 
        cat("======================\n", n1str,n2str,  "   ", nstr, "\n")
        
        res = adonis2(o ~ f, method = method, permutations = permutations)
        p = ((res$aov.tab)[["Pr(>F)"]])[1]
        if (is.na(p)) p = 1.0
        
        pmat[i,j] = round(p,roundto)
        pmat[j,i] = round(p,roundto)
        pvals[q] = round(p,roundto)
        q = q + 1
        
        if (p < 0.001) {
          smat[i,j] = "***"
          smat[j,i] = "***"
        } else if (p < 0.01) {
          smat[i,j] = "**"
          smat[j,i] = "**"
        } else if (p < 0.05) {
          smat[i,j] = "*"
          smat[j,i] = "*"
        } else if (p < 0.1) {
          smat[i,j] = "."
          smat[j,i] = "."
        } else {
          smat[i,j] = "-"
          smat[j,i] = "-"
        }
        
        print(res)
        cat("======================\n\n")
      }
    }
    for (x in 1:levs_len) {
      smat[x,x] = "-"
      pmat[x,x] = 1.0
    }
    
    qvals = p.adjust(pvals, method="fdr")
    qvals = round(qvals,roundto)
    q=1
    for (i in 1:(levs_len-1)) {
      for (j in (i + 1):levs_len) {
        pmat[i, j] = qvals[q]
        if (qvals[q] < 0.001) {
          smat[i,j] = "***"
        } else if (qvals[q] < 0.01) {
          smat[i,j] = "**"
        } else if (qvals[q] < 0.05) {
          smat[i,j] = "*"
        } else if (qvals[q] < 0.1) {
          smat[i,j] = "."
        } else {
          smat[i,j] = "-"
        }
        q = q + 1
      }
    }
    cat("")
    cat("P-value and FDR matrix:\n")
    print(pmat)
    cat("")
    cat("    P = lower left triangle\n")
    cat("    FDR = upper right triangle\n")
    cat("======================\n\n")
    cat("Symbol matrix:\n")
    print(noquote(smat))
    cat("======================\n\n")
    
  }
  
  if (filen != "")  sink()
  
  if (plotp == TRUE && levs_len > 2) otu.plot_pvalue_heat(pmat, maxp = 0.05)
  return(list(p=pmat,sym=smat))
}



otu.mul_corr <- function(frame = frame, pred = pred, method="spearman", rm.other = FALSE)
{
  if (nrow(frame) != length(pred)) stop("Outcome and predictor variables must have same lengths")
  if (rm.other == TRUE) {
    frame = frame[, colnames(frame) != "Other"]
  }
  n <- ncol(frame)
  
  if (is.numeric(pred) == FALSE) pred <- as.numeric(pred)

  col_name <- c(1:n)
  spearman_p_value <- c(1:n)
  spearman_rho <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value <- c(1:n)
  
  for (i in 1:n)
  {
    vec = frame[,i]
    sp <- cor.test(frame[,i], pred, method=method)
    spearman_p_value[i] = signif(sp$p.value, 2)
    spearman_rho[i] = sp$estimate
  }
  
  fdr_p_value <- p.adjust(spearman_p_value, method = "fdr")
  fdr_p_value <- signif(fdr_p_value, 2)
  return(list(OTU = colnames(frame), p_value = spearman_p_value, fdr = fdr_p_value, rho = spearman_rho))
}


mul_krusk <- function(frame = frame, factor = factor, verbose = v, correct = TRUE, exact = TRUE)
{
	if (nrow(frame) != length(factor)) stop("Frame and factor must have same lengths")
	
	n <- ncol(frame)
	
  
	if (is.factor(factor) == 0) factor <- as.factor(factor)
	tmp = factor(factor)
  factor = tmp
	levels <- levels(factor) #Track number of levels in factor:
	level_count <- length(levels) #Track number of levels in factor: 
	
	if (level_count < 2) stop("Factor must have > 1 levels")
	
	col_name <- c(1:n)
	krusk_p_value <- c(1:n)
	data_name <- c(1:n)
  fdr_p_value <- c(1:n)
  
  
	for (i in 1:n)
	{
		vec = frame[,i]
		x = subset(vec, factor == levels[1])
		y = subset(vec, factor == levels[2])
		ks <- kruskal.test(frame[,i] ~ factor)
		krusk_p_value[i] = signif(ks$p.value, 2)
    
		#grp1_median[i] = signif(median(x, na.rm = TRUE), 2)
    #grp1_q1[i] = signif(quantile(x, probs = 0.25, na.rm = TRUE), 2)
    #grp1_q3[i] = signif(quantile(x, probs = 0.75, na.rm = TRUE), 2)
    
		#grp2_median[i] = signif(median(y, na.rm = TRUE), 2)
    #grp2_q1[i] = signif(quantile(y, probs = 0.25, na.rm = TRUE), 2)
    #grp2_q3[i] = signif(quantile(y, probs = 0.75, na.rm = TRUE), 2)
	}
	
  fdr_p_value <- p.adjust(krusk_p_value, method = "fdr")
  fdr_p_value <- signif(fdr_p_value, 2)
	#cbind(colnames(frame), krusk_p_value, fdr_p_value, grp1_median, grp1_q1, grp1_q3, grp2_median, grp2_q1, grp2_q3)
  cbind(colnames(frame), krusk_p_value, fdr_p_value)
}

mul_wilcox <- function(frame = frame, factor = factor, min_cutoff = 0.0, sort = FALSE, verbose = FALSE, correct = TRUE, exact = TRUE)
{
	if (nrow(frame) != length(factor)) stop("Frame and factor must have same lengths")
	
	n <- ncol(frame)
	if (is.factor(factor) == 0) factor <- as.factor(factor)
	tmp = factor(factor)
  factor = tmp
	levels <- levels(factor) #Track number of levels in factor:
	level_count <- length(levels) #Track number of levels in factor: 
	
	if (level_count < 2) stop("Factor must have > 1 levels")
	col_names = colnames(frame)
	outcome <- c(1:n)
	wilcox_p_value <- c(1:n)
	krusk_p_value <- c(1:n)
	data_name <- c(1:n)
  fdr_p_value <- c(1:n)
  median_iqr1 <- c(1:n)
  median_iqr2 <- c(1:n)
  mean_sd1 <- c(1:n)
  mean_sd2 <- c(1:n)
  grp1_min <- c(1:n)
  grp2_min <- c(1:n)
  grp1_max <- c(1:n)
  grp2_max <- c(1:n)
  all_min <- c(1:n)
  all_max <- c(1:n)
  grand_mean <- c(1:n)
  #print(min_cutoff)
  
	for (i in 1:n)
	{
	  vec = frame[,i]
	  x = subset(vec, factor == levels[1])
	  y = subset(vec, factor == levels[2])
	  #all_xy = c(x,y)
	  all_xy = vec
	  grp1_min[i] = signif(min(x, na.rm = TRUE), 2)
	  grp2_min[i] = signif(min(y, na.rm = TRUE), 2)
	  grp1_max[i] = signif(max(x, na.rm = TRUE), 2)
	  grp2_max[i] = signif(max(y, na.rm = TRUE), 2)
	  all_max[i] = signif(max(c(grp1_max[i], grp2_max[i]), na.rm = TRUE), 2)
	  all_min[i] = signif(min(c(grp1_min[i], grp2_min[i]), na.rm = TRUE), 2)
	  
	  median_iqr1[i] = median_iqr_str(x)
	  median_iqr2[i] = median_iqr_str(y)
	  mean_sd1[i] = mean_sd_str(x)
	  mean_sd2[i] = mean_sd_str(y)
	  grand_mean[i] = round(mean(all_xy, na.rm = TRUE), 2)
	  
	  if (all_max[i] >= min_cutoff)
	  {
	    if (length(col_names) > 0)
	      outcome[i] = col_names[i]
	    end
		  wc <- wilcox.test(frame[,i] ~ factor, correct = correct, exact = exact)
		  ks <- kruskal.test(frame[,i] ~ factor)
		  wilcox_p_value[i] = signif(wc$p.value, 2)
		  krusk_p_value[i] = signif(ks$p.value, 2)
      
	  } else {
	    if (length(col_names) > 0)
	      outcome[i] = gsub(" ", "", paste("zExcluded-", col_names[i]))
	    end
	    wilcox_p_value[i] = NA
	    krusk_p_value[i] = NA
	  }
	}
	
  fdr_p_value <- p.adjust(wilcox_p_value, method = "fdr")
  fdr_p_value <- signif(fdr_p_value, 2)
  
	#cbind(colnames(frame), wilcox_p_value, fdr_p_value, krusk_p_value, grp1_median, grp1_q1, grp1_q3, grp2_median, grp2_q1, grp2_q3)
  
  if (verbose == TRUE)
    #data_out = data.frame(outcome = outcome, median_iqr1 = as.character(median_iqr1), mean_sd1 = as.character(mean_sd1), 
     #                 median_iqr2 = median_iqr2,  mean_sd2 = mean_sd2, wilcox_p_value = as.numeric(wilcox_p_value), fdr_p_value = as.numeric(fdr_p_value), 
     #                 grp1_max = grp1_max, grp2_max = grp2_max, all_max = all_max)
    data_out = cbind(outcome, wilcox_p_value, fdr_p_value, grand_mean, median_iqr1, mean_sd1,median_iqr2, mean_sd2, grp1_max, grp2_max, all_max)
  else
    #data_out = data.frame(outcome = outcome, median_iqr1 = median_iqr1, mean_sd1 = mean_sd1, 
     #                     median_iqr2 = median_iqr2,  mean_sd2 = mean_sd2, wilcox_p_value = as.numeric(wilcox_p_value), fdr_p_value = as.numeric(fdr_p_value))
    data_out = cbind(outcome,  wilcox_p_value, fdr_p_value, grand_mean, median_iqr1, mean_sd1,median_iqr2, mean_sd2)
  end
  
  if (sort == TRUE)
    #data_out = data_out[order(data_out$fdr_p_value, na.last = TRUE), ]
    data_out = data_out[order(fdr_p_value, na.last = TRUE), ]
  end
  
  return(data_out)
}

median_iqr_str =  function(l = data)
{
  s = summary(l)
  
  #s[3] = round(s[3],2)
  
  if (s[3] >= 1.0) 
    s[3] = round(s[3],1)
  else 
    s[3] = round(s[3],2)
  end
  
  if (s[2] >= 1.0) 
    s[2] = round(s[2],1)
  else 
    s[2] = round(s[2],2)
  end
  
  if (s[5] >= 1.0) 
    s[5] = round(s[5],1)
  else 
    s[5] = round(s[5],2)
  end
  
  res1 = paste(" (",s[2], "-",s[5],")")
  res1 = gsub(" ", "", res1)
  res2 = paste(s[3],res1)
  return(res2)
}

mean_sd_str =  function(l = data)
{
  mean = mean(l, na.rm = TRUE)
  sd = sd(l, na.rm = TRUE)
  
 # mean = round(mean,2)
  
  if (mean >= 1.0) 
    mean = round(mean, 1)
  else 
    mean = round(mean,2)
  end
  
  if (sd >= 1.0) 
    sd = round(sd,1)
  else 
    sd = round(sd,2)
  end
  
  res1 = paste(" (",sd,")")
  res1 = gsub(" ", "", res1)
  res2 = paste(mean, res1)
  return(res2)
}

mul_paired_wilcox <- function(frame1 = frame2, frame2 = frame2, min_cutoff = 0.0, sort = FALSE, verbose = FALSE, correct = TRUE, exact = TRUE)
{
	if (nrow(frame1) != nrow(frame2)) stop("Both frames must have same # of rows")
	
	n <- ncol(frame1)
  col_names = colnames(frame1)
  outcome <- c(1:n)
  wilcox_p_value <- c(1:n)
  krusk_p_value <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value <- c(1:n)
  median_iqr1 <- c(1:n)
  median_iqr2 <- c(1:n)
  mean_sd1 <- c(1:n)
  mean_sd2 <- c(1:n)
  grp1_min <- c(1:n)
  grp2_min <- c(1:n)
  grp1_max <- c(1:n)
  grp2_max <- c(1:n)
  all_min <- c(1:n)
  all_max <- c(1:n)
  grand_mean <- c(1:n)
  
  for (i in 1:n)
  {
    x = frame1[,i]
    y = frame2[,i]
    all_xy = c(x,y)
    
    grp1_min[i] = signif(min(x, na.rm = TRUE), 2)
    grp2_min[i] = signif(min(y, na.rm = TRUE), 2)
    grp1_max[i] = signif(max(x, na.rm = TRUE), 2)
    grp2_max[i] = signif(max(y, na.rm = TRUE), 2)
    all_max[i] = signif(max(c(grp1_max[i], grp2_max[i]), na.rm = TRUE), 2)
    all_min[i] = signif(min(c(grp1_min[i], grp2_min[i]), na.rm = TRUE), 2)
    
    median_iqr1[i] = median_iqr_str(x)
    median_iqr2[i] = median_iqr_str(y)
    mean_sd1[i] = mean_sd_str(x)
    mean_sd2[i] = mean_sd_str(y)
    grand_mean[i] = round(mean(all_xy, na.rm = TRUE), 2)
    
    if (all_max[i] >= min_cutoff)
    {
      if (length(col_names) > 0)
        outcome[i] = col_names[i]
      end
      #wc <- wilcox.test(frame[,i] ~ factor, correct = correct, exact = exact)
      #ks <- kruskal.test(frame[,i] ~ factor)
      #wilcox_p_value[i] = signif(wc$p.value, 2)
      #krusk_p_value[i] = signif(ks$p.value, 2)
      
      wc <- wilcox.test(frame1[,i], frame2[,i], paired = TRUE, correct = correct, exact = exact)
      #ks <- kruskal.test(frame1[,i], frame2[,i], paired = TRUE)
      wilcox_p_value[i] = signif(wc$p.value, 2)
      krusk_p_value[i] = NA #signif(ks$p.value, 2)
      
    } else {
      if (length(col_names) > 0)
        outcome[i] = gsub(" ", "", paste("zExcluded-", col_names[i]))
      end
      wilcox_p_value[i] = NA
      krusk_p_value[i] = NA
    }
  }
  fdr_p_value <- p.adjust(wilcox_p_value, method = "fdr")
  fdr_p_value <- signif(fdr_p_value, 2)
 
	if (verbose == TRUE)
	  data_out = cbind(outcome, grand_mean, median_iqr1, mean_sd1,median_iqr2, mean_sd2, wilcox_p_value, fdr_p_value, grp1_max, grp2_max, all_max)
	else
	  data_out = cbind(outcome, grand_mean, median_iqr1, mean_sd1,median_iqr2, mean_sd2, wilcox_p_value, fdr_p_value)
	end
	
	if (sort == TRUE)
	  data_out = data_out[order(fdr_p_value, na.last = TRUE), ]
	end
	
	return(data_out)
}

mul_aov <- function(frame = frame, factor1 = factor1, factor2 = "")
{
  if (nrow(frame) != length(factor1)) stop("Frame and factor must have same lengths")
  
  n <- ncol(frame)

  if (is.factor(factor1) == FALSE) factor1 <- as.factor(factor1)
  #tmp = factor(factor)
  #factor = tmp
  levels <- levels(factor1) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 
  #print(level_count)
  if (level_count < 2) stop("Factor must have > 1 levels")
  
  col_name <- c(1:n)
  p_value1 <- c(1:n)
  p_value2 <- c(1:n)
  data_name <- c(1:n)
  fdr_value1 <- c(1:n)
  fdr_value2 <- c(1:n)
  
  for (i in 1:n)
  {
    #vec = frame[,i]
    #x = subset(vec, factor == levels[1])
    #y = subset(vec, factor == levels[2])
    
    if (factor2 == "") {
      result = summary(aov(frame[,i] ~ factor1, na.rm = TRUE))
    } else {
      result = summary(aov(frame[,i] ~ factor1 + factor2, na.rm = TRUE))
    }

    p1 = result[[1]][["Pr(>F)"]][[1]]
    p_value1[i] = signif(p1, 2)
    
    if (factor2 != "") {
      p2 = result[[1]][["Pr(>F)"]][[2]]
      p_value2[i] = signif(p2, 2)
    }
      
    #grp1_median[i] = signif(median(x, na.rm = TRUE), 2)
    #grp1_q1[i] = signif(quantile(x, probs = 0.25, na.rm = TRUE), 2)
    #grp1_q3[i] = signif(quantile(x, probs = 0.75, na.rm = TRUE), 2)
    
    #grp2_median[i] = signif(median(y, na.rm = TRUE), 2)
    #grp2_q1[i] = signif(quantile(y, probs = 0.25, na.rm = TRUE), 2)
    #grp2_q3[i] = signif(quantile(y, probs = 0.75, na.rm = TRUE), 2)
  }
  
  fdr_value1 <- p.adjust(p_value1, method = "fdr")
  fdr_value1 <- signif(fdr_value1, 2)
  #cbind(colnames(frame), krusk_p_value, fdr_p_value, grp1_median, grp1_q1, grp1_q3, grp2_median, grp2_q1, grp2_q3)
  
  if (factor2 == "") {
    r =cbind(colnames(frame), p_value1, fdr_value1)
  } else {
    fdr_value2 <- p.adjust(p_value2, method = "fdr")
    fdr_value2 <- signif(fdr_value2, 2)
    
    r =cbind(colnames(frame), p_value1, fdr_value1, p_value2, fdr_value2)
  }
  
  return(r)
}

mul_aov_tukey <- function(frame = frame, factor = factor)
{
  factor = droplevels(factor)  # Exclude any unused levels
  if (nrow(frame) != length(factor)) stop("Frame and factor must have same lengths")
  
  n <- ncol(frame)
  
  if (is.factor(factor) == FALSE) factor <- as.factor(factor)

  levels <- levels(factor) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 

  if (level_count < 2) stop("Factor must have > 1 levels")
  
  cols = (level_count * (level_count - 1))/2
  tuk_mat = matrix(, nrow = n, ncol = cols)
  
  p_value <- c(1:n)
  fdr_value <- c(1:n)
  
  for (i in 1:n) {
    
    aov_fit = aov(frame[,i] ~ factor, na.rm = TRUE)
    aov_sum = summary(aov_fit)
    p = aov_sum[[1]][["Pr(>F)"]][[1]]
    p_value[i] = signif(p, 2)
    
    if (level_count > 2) {
      tuk_res_mat = (TukeyHSD(aov_fit))[[1]]
      tuk_mat[i,] = signif(tuk_res_mat[,4], 2)
      colnames(tuk_mat) = rownames(tuk_res_mat)
    }
  }
  
  #print(tuk_mat)
  
  fdr_value <- p.adjust(p_value, method = "fdr")
  fdr_value <- signif(fdr_value, 2)
  
  if (level_count > 2) {
    r = cbind(colnames(frame), p_value, fdr_value, tuk_mat)
  } else {
    r = cbind(colnames(frame), p_value, fdr_value)
  }
  return(r)
}

old_mul_lm <- function(frame = frame, factor1 = factor1, factor2 = "")
{
  if (nrow(frame) != length(factor1)) stop("Frame and factor must have same lengths")
  
  n <- ncol(frame)
  if (is.factor(factor1) == FALSE) factor1 <- as.factor(factor1)

  levels <- levels(factor1) #Track number of levels in factor:
  level_count <- length(levels) #Track number of levels in factor: 

  if (level_count < 2) stop("Factor must have > 1 levels")
  
  col_name <- c(1:n)
  lm_p_value1 <- c(1:n)
  lm_p_value2 <- c(1:n)
  data_name <- c(1:n)
  fdr_p_value1 <- c(1:n)
  fdr_p_value2 <- c(1:n)
  
  for (i in 1:n)
  {
    #vec = frame[,i]
    #x = subset(vec, factor == levels[1])
    #y = subset(vec, factor == levels[2])
    
    if (length(factor2) <= 1) {
      result = coefficients(summary(lm(frame[,i] ~ factor1, na.action = na.omit)))
    } else {
      result = coefficients(summary(lm(frame[,i] ~ factor1 + factor2, na.action = na.omit)))
    }
   # print(result)
    lm_p_value1[i] = signif(result[2,4], 2)
    
    #grp1_median[i] = signif(median(x, na.rm = TRUE), 2)
    #grp1_q1[i] = signif(quantile(x, probs = 0.25, na.rm = TRUE), 2)
    #grp1_q3[i] = signif(quantile(x, probs = 0.75, na.rm = TRUE), 2)
    
    #grp2_median[i] = signif(median(y, na.rm = TRUE), 2)
    #grp2_q1[i] = signif(quantile(y, probs = 0.25, na.rm = TRUE), 2)
    #grp2_q3[i] = signif(quantile(y, probs = 0.75, na.rm = TRUE), 2)
  }
  
  fdr_p_value1 <- p.adjust(lm_p_value1, method = "fdr")
  fdr_p_value1 <- signif(fdr_p_value1, 2)
  #cbind(colnames(frame), krusk_p_value, fdr_p_value, grp1_median, grp1_q1, grp1_q3, grp2_median, grp2_q1, grp2_q3)
  cbind(colnames(frame), lm_p_value1, fdr_p_value1)
}

median_iqr_str <- function(l = data, rnd = -1)
{
  qs = quantile(l, probs=c(0.25, 0.50, 0.75), na.rm = TRUE)
  if (rnd != -1) {
    qs = round(qs, rnd)
    res = paste(format(qs[2], nsmall=rnd)," (", format(qs[1], nsmall=rnd), "-", format(qs[3], nsmall=rnd), ")", sep="")
  } else {
    res = paste(format(qs[2])," (", format(qs[1]), "-", format(qs[3]), ")", sep="")
  }
  return(res)
}

mean_sd_str <- function(l = data, rnd = -1)
{
  mn = mean(l, na.rm = TRUE)
  sd = sd(l,  na.rm = TRUE)
  if (rnd != -1) {
    mn = round(mn, rnd)
    sd = round(sd, rnd)
    res = paste(format(mn, nsmall=rnd)," (", format(sd, nsmall=rnd), ")", sep="")
  } else {
    res = paste(format(mn)," (", format(sd), ")", sep="")
  }
  
  return(res)
}

oldmean_sd_str =  function(l = data)
{
  mean = mean(l, na.rm = TRUE)
  sd = sd(l, na.rm = TRUE)
  
  # mean = round(mean,2)
  
  if (mean >= 1.0) 
    mean = round(mean, 1)
  else 
    mean = round(mean,2)
  end
  
  if (sd >= 1.0) 
    sd = round(sd,1)
  else 
    sd = round(sd,2)
  end
  
  res1 = paste(" (",sd,")")
  res1 = gsub(" ", "", res1)
  res2 = paste(mean, res1)
  return(res2)
}