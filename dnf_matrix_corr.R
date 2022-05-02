#************************
#
#
# v.0.1.0  
# daniel.frank@ucdenver.edu
# 24 April 2017
#
#
#************************
library("car")  #For Type II anova p-values

get_corr_mat <- function(mat1 = matrix, 
                         mat2 = matrix, 
                         method = c("pearson", "kendall", "spearman"),
                         maxp = 0.05,
                         diag = TRUE)
  
{
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  est_mat = matrix(nrow = ncol1, ncol = ncol2)
  p_mat = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat = matrix(nrow = ncol1, ncol = ncol2)

  colnames(est_mat) = colnames(mat2)
  rownames(est_mat) = colnames(mat1)
  colnames(p_mat) = colnames(mat2)
  rownames(p_mat) = colnames(mat1)
  colnames(logp_mat) = colnames(mat2)
  rownames(logp_mat) = colnames(mat1)
  colnames(sym_mat) = colnames(mat2)
  rownames(sym_mat) = colnames(mat1)

  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      nas = !is.na(mat1[,i]) & !is.na(mat2[,j])
      if (((sum(!is.na(mat2[,j])) < 1) | (sum(!is.na(mat1[,i])) < 1)) | (sum(nas) <= 2)) {
        # est_mat[i,j] = 0
        est_mat[i,j] = NA
        p_mat[i,j] = 1
        logp_mat[i,j] = 0
        sym_mat[i,j] = ""
        next
      }
      
      if ((diag == FALSE) & (i == j)) {
        est_mat[i,j] = 0
        p_mat[i,j] = 1
        logp_mat[i,j] = 0
        sym_mat[i,j] = ""
      } else {
        result = cor.test(mat1[,i], mat2[,j], method = method)
        
        est_mat[i, j] = result$estimate
        p_mat[i, j] = result$p.value
        lp1 = get_logp(result$p.value, maxp)
        logp_mat[i, j] =lp1
        sym_mat[i,j] = get_sym(lp1)
      }
     }
  }
  return(list(est = est_mat, p = p_mat, logp = logp_mat, sym = sym_mat))
}

get_corr_edgelist <- function(mat1 = matrix, 
                         mat2 = matrix,
                         method = c("pearson", "kendall", "spearman"), 
                         rm.other = FALSE)
  
{
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  totrows = ncol1 * ncol2
  #cat("totrows:",totrows, "\n")
  index1 = vector(length = totrows)
  #cat("len index1 before:",length(index1), "\n")
  index2 = vector(length = totrows)
  estimate = vector(length = totrows)
  pvalue = vector(length = totrows)
  var1 = vector(length = totrows)
  var2 = vector(length = totrows)
  tot = 0
  row = 1
  for (i in 1:ncol1) {
    cat("i:",i,"of", ncol1,"\n")
    for (j in 1:ncol2) {
      
      nas = !is.na(mat1[,i]) & !is.na(mat2[,j])
      #print(mat1[,i])
      #print(mat2[,j])
      #print(sum(nas))
      if (((sum(!is.na(mat2[,j])) < 1) | (sum(!is.na(mat1[,i])) < 1)) | (sum(nas) <= 2)) {
         next
      }
      #if (i > 102) cat("j:",j,"of", ncol2,"\n")
      
      result = cor.test(mat1[,i], mat2[,j], method = method)
    
      index1[row] = i
      index2[row] =  j
      estimate[row] = result$estimate
      if (is.na(result$p.value)) {
        pvalue[row] = 1.0
      } else {
        pvalue[row] = result$p.value
      }
      var1[row] = (colnames(mat1))[i]
      var2[row] = (colnames(mat2))[j]
      row = row + 1
    }
  }
  #cat("length index1:", length(index1), "\n")
  res_df = data.frame(index1, index2, estimate, pvalue, var1, var2)
  #cat("total:", row, "\n")
  return(res_df)
}

get_lm_edgelist <- function(mat1 = matrix, 
                              mat2 = matrix,
                              factor = "",
                              showinteraction = FALSE,
                              rm.other = FALSE)
  
{
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  totrows = ncol1 * ncol2
  index1 = vector(length = totrows)
  index2 = vector(length = totrows)
  estimate = vector(length = totrows)
  pvalue = vector(length = totrows)
  var1 = vector(length = totrows)
  var2 = vector(length = totrows)
  
  row = 1
  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      if (length(factor) <= 1) {
        result = coefficients(summary(lm(mat1[,i] ~ mat2[,j])))
        estimate[row] = result[2,3]
        if (is.na(result[2,4])) {
          pvalue[row] = 1.0
        } else {
          pvalue[row] = result[2,4]
        }
        
      } else {
        #result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] * factor)))
        #if (showinteraction == FALSE) {
        #  estimate[row] = result[2,3]
        #  if (is.na(result[2,4])) {
        #    pvalue[row] = 1.0
        #  } else {
        #    pvalue[row] = result[2,4]
        #  }
        #} else if (showinteraction == TRUE){
        #  estimate[row] = result[4,3]
        #  if (is.na(result[4,4])) {
        #    pvalue[row] = 1.0
        #  } else {
        #    pvalue[row] = result[4,4]
        #  }
        #}
        result = coefficients(summary(lm(mat1[,i] ~ factor * mat2[,j])))
        
        if (showinteraction == FALSE) {
          estimate[row] = result[3,3]
          if (is.na(result[3,4])) {
            pvalue[row] = 1.0
          } else {
            pvalue[row] = result[3,4]
          }
        } else if (showinteraction == TRUE){
          estimate[row] = result[4,3]
          if (is.na(result[4,4])) {
            pvalue[row] = 1.0
          } else {
            pvalue[row] = result[4,4]
          }
        }
      }
      
      index1[row] = i
     index2[row] =  j
      #estimate[row] = result$estimate
     # if (is.na(result$p.value)) {
     #   pvalue[row] = 1.0
     # } else {
     #   pvalue[row] = result$p.value
     # }
      var1[row] = (colnames(mat1))[i]
      var2[row] = (colnames(mat2))[j]
      row = row + 1
    }
  }
  res_df = data.frame(index1, index2, estimate, pvalue, var1, var2)
  
  return(res_df)
}

matrix_corr =  function(mat1 = matrix, 
                         mat2 = matrix, 
                         method = c("pearson", "kendall", "spearman"), 
                         heatmap = "TRUE", 
                         cluster = c("row", "col", "column", "both", "none"), 
                         pvalue = "FALSE", 
                         maxp = 0.1,
                         diag = "TRUE",
                         filen = "",
                         main = "",
                         pdfx = 10, pdfy = 10, rev_col = FALSE, srtCol = 90, keysize = 1, 
                         plotsymbols = TRUE, rm.other = FALSE, notecex = 1, ColSideColors = NULL,
                         flt_col_maxp = 1.0, flt_row_maxp = 1.0)
{
  
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  corr_res = get_corr_mat(mat1, mat2, method = method,
                          maxp = maxp, diag = diag)

  corr_res = shrink_corr_mat(corr_res, flt_col_maxp = flt_col_maxp, flt_row_maxp = flt_row_maxp)
  
  out_mat = corr_res$est
  if (pvalue == TRUE) {
    cluster = "none"
    out_mat = corr_res$p
  }
  
  if (heatmap == "TRUE") {
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
    }
    
    if (pvalue == "TRUE") {
      if (rev_col == TRUE) {
        colpan=colorpanel(100,"white","blue","red")
      } else {
        colpan=colorpanel(100,"white","blue","darkblue")
      }
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat,col=colpan,trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=corr_res$sym, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = TRUE)
      } else {
        heatmap.2(out_mat,col=colpan,trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=corr_res$sym, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = TRUE, ColSideColors = ColSideColors)
      }
    } else {
      keytitle = ""
      if (method == "pearson") {
        keytitle = "Pearson R"
      } else if (method == "spearman") {
        keytitle = "Spearman rho"
      }
      if (rev_col == TRUE) {
        colpan=colorpanel(100,"red","white","blue")
      } else {
        colpan=colorpanel(100,"blue","white","red")
      }
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat,col=colpan, trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=corr_res$sym, notecol = "white", keysize=keysize, key.xlab = keytitle, 
                  key.title = NA, notecex = notecex, symkey = TRUE, symbreaks = TRUE)
      } else {
        heatmap.2(out_mat,col=colpan, trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=corr_res$sym, notecol = "white", keysize=keysize, key.xlab = keytitle, 
                  key.title = NA, notecex = notecex, symkey = TRUE, symbreaks = TRUE, ColSideColors = ColSideColors)
      }
    }
    
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
  }
 
  res = list(est = corr_res$est, logp = corr_res$logp, sym = corr_res$sym)
  return(res)
}

shrink_corr_mat =  function(corr_res, flt_col_maxp = 1.0, flt_row_maxp = 1.0)
{
  out_mat1 = corr_res$est
  #print(out_mat1)
  cat("Initial NCOL: ", ncol(out_mat1), "\n")
  cat("Initial NROW: ", nrow(out_mat1), "\n")
  
  logp_mat1 = corr_res$logp
  sym_mat1 = corr_res$sym
  
  rows = nrow(out_mat1)
  inc_rows = rep(TRUE, rows)
 
  cols = ncol(out_mat1)
  inc_cols = rep(TRUE, cols)
 
  log_col_maxp = -1 * log10(flt_col_maxp)
  log_row_maxp = -1 * log10(flt_row_maxp)
  #print(log_col_maxp)
  #print(log_row_maxp)
  
  for (i in 1:rows) {
    max1 = max(logp_mat1[i, ])
    if (max1 < log_row_maxp) {
      inc_rows[i] = FALSE
    }
    if (sum(!is.na(out_mat1[i, ])) < 1) inc_rows[i] = FALSE
    
  }
  #  print( sum(inc_rows))
  
  #Only alter matrices if at least two rows will be retained
  if (sum(inc_rows) > 1) {
    out_mat1 = out_mat1[inc_rows == TRUE, ]
    logp_mat1 = logp_mat1[inc_rows == TRUE, ]
    sym_mat1 = sym_mat1[inc_rows == TRUE, ]    
  }
  
  for (i in 1:cols) {
    max1 = max(logp_mat1[, i])
    if (max1 < log_col_maxp) inc_cols[i] = FALSE
    if (sum(!is.na(out_mat1[ ,i])) < 1) inc_cols[i] = FALSE
  }
  
  #Only alter matrices if at least two columns will be retained
  if (sum(inc_cols) > 1) {
    out_mat1 = out_mat1[, inc_cols == TRUE]
    logp_mat1 = logp_mat1[, inc_cols == TRUE]
    sym_mat1 = sym_mat1[, inc_cols == TRUE]
  }
  
  cat("Final NCOL: ", ncol(out_mat1), "\n")
  cat("Final NROW: ", nrow(out_mat1), "\n")
  
  res = list(est = out_mat1, logp = logp_mat1, sym = sym_mat1)
  return(res)
}


matrix_corr0 =  function(mat1 = matrix, 
                        mat2 = matrix, 
                        method = c("pearson", "kendall", "spearman"), 
                        heatmap = "FALSE", 
                        cluster = c("row", "col", "column", "both", "none"), 
                        pvalue = "FALSE", 
                        maxp = 0.05,
                        diag = "TRUE",
                        filen = "",
                        main = "",
                        pdfx = 10, pdfy = 10, rev_col = FALSE, srtCol = 90, keysize = 1, 
                        plotsymbols = TRUE, rm.other = FALSE, notecex = 1, ColSideColors = NULL)
{
  
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  sym_mat = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat = matrix(nrow = ncol1, ncol = ncol2)
  res = get_corr_mat(mat1, mat2, method = method)
  est_mat = res$est
  p_mat = res$p
  
  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      if ((diag == FALSE) & (i == j)) {
        est_mat[i,j] = NA
        p_mat[i,j] = NA
        logp_mat[i,j] = NA
      } else {
        #Set min low value so -log(pvalue) range isn't so great
        if (p_mat[i, j] < 0.000001) {
          p_mat[i, j] = 0.000001
        } 
      
        logp_mat[i, j] = -1*log10(p_mat[i, j])
        
        if (logp_mat[i,j] < -1*log10(maxp)) logp_mat[i,j] = 0
        
        if (logp_mat[i,j] >= -log10(0.001) & plotsymbols == TRUE) {
          sym_mat[i, j] = '#'
        } else if (logp_mat[i,j] >= -log10(0.01) & plotsymbols == TRUE) {
          sym_mat[i, j] = '+'
        } else if (logp_mat[i,j] >= -log10(0.05) & plotsymbols == TRUE) {
          #sym_mat[i, j] = '+'
          sym_mat[i, j] = '*'
        }
      }
    }
  }
 
  out_mat = est_mat
  if (pvalue == TRUE) {
    cluster = "none"
    out_mat = p_mat
  }
  
  if (heatmap == "TRUE") {
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
    }
 
    if (pvalue == "TRUE") {
      if (rev_col == TRUE) {
        colpan=colorpanel(100,"white","blue","red")
      } else {
        colpan=colorpanel(100,"white","blue","darkblue")
      }
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat,col=colpan,trace="none", Colv=colv, Rowv=rowv, 
              main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
              key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = TRUE)
      } else {
        heatmap.2(out_mat,col=colpan,trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = TRUE, ColSideColors = ColSideColors)
      }
    } else {
      keytitle = ""
      if (method == "pearson") {
        keytitle = "Pearson R"
      } else if (method == "spearman") {
        keytitle = "Spearman rho"
      }
      if (rev_col == TRUE) {
        colpan=colorpanel(100,"red","white","blue")
      } else {
        colpan=colorpanel(100,"blue","white","red")
      }
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat,col=colpan, trace="none", Colv=colv, Rowv=rowv, 
              main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab = keytitle, 
              key.title = NA, notecex = notecex, symkey = TRUE, symbreaks = TRUE)
      } else {
        heatmap.2(out_mat,col=colpan, trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab = keytitle, 
                  key.title = NA, notecex = notecex, symkey = TRUE, symbreaks = TRUE, ColSideColors = ColSideColors)
      }
    }
    
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
  }
  return(out_mat)
}

dualmatrix_corr =  function(mat1 = matrix, 
                        mat2 = matrix, 
                        mat3 = matrix,
                        mat4 = matrix,
                        method = c("pearson", "kendall", "spearman"), 
                        heatmap = "FALSE", 
                        cluster = c("row", "col", "column", "both", "none"), 
                        pvalue = "FALSE", 
                        maxp = 0.05,
                        diag = "TRUE",
                        main = "",
                        filen = "",
                        pdfx = 10, pdfy = 10, rev_col = FALSE,
                        srtCol = 90, keysize = 1, plotsymbols = TRUE, rm.other = FALSE, notecex = 1)
{
  if (cluster == "col" || cluster == "column" || cluster == "both") {
    print("Clustering only by row to separate matrices.")
  }
  if ((nrow(mat1) != nrow(mat2)) || (nrow(mat3) != nrow(mat4))) stop("All matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
    mat3 = mat3[, colnames(mat1) != "Other"]
    mat4 = mat4[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  ncol3 = ncol(mat3)
  ncol4 = ncol(mat4)

  out_mat = matrix(nrow = ncol1, ncol = ncol2 + ncol4 + 1)
  p_mat = matrix(nrow = ncol1, ncol = ncol2 + ncol4 + 1)
  sym_mat = matrix(nrow = ncol1, ncol = ncol2 + ncol4 + 1)
  col_names = c(colnames(mat2)," ",colnames(mat4))
  #print(col_names)
  colnames(out_mat) = col_names
  rownames(out_mat) = colnames(mat1)
  colnames(p_mat) = col_names
  rownames(p_mat) = colnames(mat1)
  
  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      if ((diag == "TRUE") || (i != j)) {
        result = cor.test(mat1[,i], mat2[,j], method = method)
        out_mat[i, j] = result$estimate

        p = result$p.value
        if (is.na(p) == TRUE) {
          p = 1.0
        }
        else if (p < 0.0000001) {
          p = 0.000001
        } 
        lp = -1*log10(p)
        if (lp < -1*log10(maxp)) lp = 0
        p_mat[i, j] = lp
        
        if (lp >= -1*log10(0.01) & plotsymbols == TRUE) {
          sym_mat[i, j] = '#'
        } 
        else if (lp >= -log10(0.05) & plotsymbols == TRUE) {
          sym_mat[i, j] = '+'
        }
      }
    }
  }
  
  
  for (i in 1:ncol3) {
    for (j in 1:ncol4) {
      if ((diag == "TRUE") || (i != j)) {
        result = cor.test(mat3[,i], mat4[,j], method = method)
        out_mat[i, j+ncol2+1] = result$estimate

        p = result$p.value
          
        if (is.na(p) == TRUE) {
           p = 1.0
        }
        else if (p < 0.0000001) {
          p = 0.000001
        } 
        lp = -1*log10(p)
        if (lp < -1*log10(maxp)) lp = 0
        p_mat[i, j+ncol2+1] = lp
          
        if (lp >= -log10(0.01) & plotsymbols == TRUE) {
          sym_mat[i, j+ncol2+1] = '#'
        } else if (lp >= -log10(0.05) & plotsymbols == TRUE) {
          sym_mat[i, j+ncol2+1] = '+'
        }
      }
    }
  }
  
  if (pvalue == "TRUE") {
    cluster = "none"
    out_mat = p_mat
  }
  
  
  if (heatmap == "TRUE") {
    m = out_mat
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
      #print("both")
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
      #print("row")
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
      #print("col")
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
      #print("row")
    }
    
    #if (pvalue == "TRUE") {
    #  if (rev_col == FALSE) {
    #    heatmap.2(m,col=colorpanel(100,"white","blue","red"),trace="none", Colv=colv, Rowv=rowv, main=main)
    #  } else {
    #    heatmap.2(m,col=colorpanel(100,"white","blue","darkblue"),trace="none", Colv=colv, Rowv=rowv, main=main)
    #  }
   # } else {
    #  if (rev_col == FALSE) {
    #    heatmap.2(m,col=colorpanel(100,"blue","white", "red"),trace="none", Colv=colv, Rowv=rowv, main=main)
    #  } else {
    #    heatmap.2(m,col=colorpanel(100,"red","white", "blue"),trace="none", Colv=colv, Rowv=rowv, main=main)
    #  }
   # }
  
    
    if (pvalue == "TRUE") {
      if (rev_col == FALSE) {
        heatmap.2(m,col=colorpanel(100,"white","blue","red"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE)
      } else {
        heatmap.2(m,col=colorpanel(100,"white","blue","darkblue"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE)
      }
    } else {
      keytitle = ""
      if (method == "pearson") {
        keytitle = "Pearson R"
      } else if (method == "spearman") {
        keytitle = "Spearman rho"
      }
    
      if (rev_col == FALSE) {
        heatmap.2(m,col=colorpanel(100,"blue","white", "red"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab = keytitle, 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE)
      } else {
        heatmap.2(m,col=colorpanel(100,"red","white", "blue"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab = keytitle, 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE)
      }
    }
    
    #heatmap.2(m,trace="none", Colv=colv, Rowv=rowv)
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
  }
 
  return(out_mat)
}


trimatrix_corr =  function(mat1 = matrix, 
                        mat2 = matrix, 
                        mat3 = matrix,
                        mat4 = matrix,
                        mat5 = matrix,
                        mat6 = matrix,
                        method = c("pearson", "kendall", "spearman"), 
                        heatmap = "FALSE", 
                        cluster = c("row", "col", "column", "both", "none"), 
                        pvalue = "FALSE", 
                        maxp = 0.05,
                        diag = "TRUE",
                        lower = "FALSE",
                        filen = "",
                        pdfx = 8, pdfy = 30, rev_col = FALSE)
{
  if (cluster == "col" || cluster == "column" || cluster == "both") {
    print("Clustering only by row to separate matrices.")
  }
  if ((nrow(mat1) != nrow(mat2)) || (nrow(mat3) != nrow(mat4))) stop("All matrices must have the same number of rows")
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  ncol3 = ncol(mat3)
  ncol4 = ncol(mat4)
  ncol5 = ncol(mat5)
  ncol6 = ncol(mat6)
  
  out_mat = matrix(nrow = ncol1, ncol = ncol2 + ncol4 + ncol6 + 2)
  col_names = c(colnames(mat2)," ",colnames(mat4), " ",colnames(mat6))
  #print(col_names)
  colnames(out_mat) = col_names
  rownames(out_mat) = colnames(mat1)
  
  for (i in 1:ncol1) {
    k= ncol2
    if (lower == "TRUE") k = i
    for (j in 1:k) {
      if ((diag == "TRUE") || (i != j)) {
        result = cor.test(mat1[,i], mat2[,j], method = method)
        out_mat[i, j] = result$estimate
        if (pvalue == "TRUE") {
          cluster = "none"
          p = result$p.value
        
          if (is.na(p) == TRUE) {
            p = 1.0
          }
          else if (p < 0.0000001) {
            p = 0.000001
          } 
          lp = -1*log10(p)
          if (lp < -1*log10(maxp)) lp = NA
          out_mat[i, j] = lp
        }
      }
    }
  }
  
  for (i in 1:ncol3) {
    k= ncol4
    if (lower == "TRUE") k = i
    
    for (j in 1:k) {
      if ((diag == "TRUE") || (i != j)) {
        result = cor.test(mat3[,i], mat4[,j], method = method)
        out_mat[i, j+ncol2+1] = result$estimate
        if (pvalue == "TRUE") {
          cluster = "none"
          p = result$p.value
          
          if (is.na(p) == TRUE) {
            p = 1.0
          }
          else if (p < 0.0000001) {
            p = 0.000001
          } 
          lp = -1*log10(p)
          if (lp < -1*log10(maxp)) lp = NA
          out_mat[i, j+ncol2+1] = lp
        }
      }
    }
  }
  
  for (i in 1:ncol5) {
    k= ncol6
    if (lower == "TRUE") k = i
    
    for (j in 1:k) {
      
      if ((diag == "TRUE") || (i != j)) {
        result = cor.test(mat5[,i], mat6[,j], method = method)
        out_mat[i, j+ ncol2 + ncol4+2] = result$estimate
        if (pvalue == "TRUE") {
          cluster = "none"
          p = result$p.value
          
          if (is.na(p) == TRUE) {
            p = 1.0
          }
          else if (p < 0.0000001) {
            p = 0.000001
          } 
          lp = -1*log10(p)
          if (lp < -1*log10(maxp)) lp = NA
          out_mat[i, j+ ncol2 + ncol4+2] = lp
        }
      }
    }
  }
  
  if (heatmap == "TRUE") {
    m = out_mat
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
      #print("both")
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
      #print("row")
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
      #print("col")
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
      #print("row")
    }
    if (pvalue == "TRUE") {
      heatmap.2(m,col=colorpanel(100,"white","blue","darkblue"),trace="none", Colv=colv, Rowv=rowv)
      #heatmap.2(m,col=cm.colors,trace="none", Colv=colv, Rowv=rowv)
    } else {
      if (rev_col == FALSE) {
        heatmap.2(m,col=colorpanel(100,"blue","white", "red"),trace="none", Colv=colv, Rowv=rowv, keysize=1)
      } else {
        heatmap.2(m,col=colorpanel(100,"red","white", "blue"),trace="none", Colv=colv, Rowv=rowv)
      }
    }
  
    #heatmap.2(m,trace="none", Colv=colv, Rowv=rowv)
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
  }
  
  return(out_mat)
}

matrix_lm0 =  function(mat1 = matrix, 
                      mat2 = matrix, 
                      factor = "",
                      #covar = "",
                      heatmap = "FALSE", 
                      cluster = c("row", "col", "column", "both", "none"),
                      interaction = "FALSE",
                      pvalue = "FALSE", 
                      maxp = 0.05,
                      diag = "TRUE",
                      filen = "",
                      main = "",
                      pdfx = 10, pdfy = 10, rev_col = FALSE, srtCol = 90, keysize = 1, 
                      plotsymbols = TRUE, rm.other = FALSE, notecex = 1, ColSideColors = NULL , RowSideColors = NULL )
{
  
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  out_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  
  out_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  
  out_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  
  colnames(out_mat1) = colnames(mat2)
  rownames(out_mat1) = colnames(mat1)
  colnames(p_mat1) = colnames(mat2)
  rownames(p_mat1) = colnames(mat1)
  colnames(logp_mat1) = colnames(mat2)
  rownames(logp_mat1) = colnames(mat1)
  
  colnames(out_mat2) = colnames(mat2)
  rownames(out_mat2) = colnames(mat1)
  colnames(p_mat2) = colnames(mat2)
  rownames(p_mat2) = colnames(mat1)
  colnames(logp_mat2) = colnames(mat2)
  rownames(logp_mat2) = colnames(mat1)
  
  colnames(out_mat3) = colnames(mat2)
  rownames(out_mat3) = colnames(mat1)
  colnames(p_mat3) = colnames(mat2)
  rownames(p_mat3) = colnames(mat1)
  colnames(logp_mat3) = colnames(mat2)
  rownames(logp_mat3) = colnames(mat1)
  
  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      print(j)
      if ((diag == "TRUE") || (i != j)) {
        
        if (length(factor) <= 1) {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j])))
          
          out_mat1[i, j] = result[2,3]
          p_mat1[i, j] = result[2,4]
          lp1 = get_logp(result[2,4], maxp)
          logp_mat1[i, j] =lp1
          
          sym_mat1[i,j] = get_sym(lp1)
          
        } else if (interaction == FALSE) {
          
          #print(mat1[,i])
          #print(mat2[,j])
          #print(factor)
          #print("\n")
          
          factest = factor[!is.na(mat2[,j])]
          print(factest)
          print(levels(factest))
          print(length(levels(droplevels(factest))))
          print("\n")
          
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] + factor)))
          
          out_mat1[i, j] = result[2,3]
          p_mat1[i, j] = result[2,4]
          lp1 = get_logp(result[2,4], maxp)
          logp_mat1[i, j] =lp1
          
          if (nrow(result) >= 3) {
            out_mat2[i, j] = result[3,3]
            p_mat2[i, j] = result[3,4]
            lp2 = get_logp(result[3,4], maxp)
            logp_mat2[i, j] =lp2
          } else {
            out_mat2[i, j] = 0
            p_mat2[i, j] = NA
            lp2 = get_logp(NA, maxp)
            logp_mat2[i, j] =lp2
            
          }
          sym_mat1[i,j] = get_sym(lp1)
          sym_mat2[i,j] = get_sym(lp2)
        } else {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] * factor)))
         
          out_mat1[i, j] = result[2,3]
          p_mat1[i, j] = result[2,4]
          lp1 = get_logp(result[2,4], maxp)
          logp_mat1[i, j] =lp1
          
          out_mat2[i, j] = result[3,3]
          p_mat2[i, j] = result[3,4]
          lp2 = get_logp(result[3,4], maxp)
          logp_mat2[i, j] =lp2
          
          out_mat3[i, j] = result[4,3]
          p_mat3[i, j] = result[4,4]
          lp3 = get_logp(result[4,4], maxp)
          logp_mat3[i, j] =lp3
          
          sym_mat1[i,j] = get_sym(lp1)
          sym_mat2[i,j] = get_sym(lp2)
          sym_mat3[i,j] = get_sym(lp3)		
        }
      }
    }
  }
  
  if (pvalue == "TRUE") {
    cluster = "none"
    out_mat1 = logp_mat1
    out_mat2 = logp_mat2
    out_mat3 = logp_mat3
  }
  
  if (heatmap == "TRUE") {
   # m = out_mat1
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
      #print("both")
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
      #print("row")
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
      #print("col")
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
      #print("row")
    }
    
    if (pvalue == "TRUE") {
      if (rev_col == FALSE) {
        cols = colorpanel(100,"white","blue","red")
      } else {
        cols = colorpanel(100,"white","blue","darkblue")
      }
      xlab="-log10(P)"
    } else {
      if (rev_col == FALSE) {
        cols = colorpanel(100,"blue","white", "red")
      } else {
        cols = col=colorpanel(100,"red","white", "blue")
      }
      xlab="LM t-value"
    }
    
    symkey = TRUE
    
    #print(out_mat1)
    if (length(factor) > 1) {
      main1 = paste(main, "(Factor1)")
    } else {
      main1 = main
    }
    if (is.null(ColSideColors)) {
      heatmap.2(out_mat1,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                main=main1, srtCol = srtCol, density.info = "none", cellnote=sym_mat1, notecol = "white", keysize=keysize, key.xlab=xlab, 
                key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey)
    } else {
      heatmap.2(out_mat1,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                main=main1, srtCol = srtCol, density.info = "none", cellnote=sym_mat1, notecol = "white", keysize=keysize, key.xlab=xlab, 
                key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors)
    }
    
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
    
    if (length(factor) > 1) {
      main2 = paste(main, "(Factor2)")
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat2,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main2, srtCol = srtCol, density.info = "none", cellnote=sym_mat2, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey)
      } else {
        heatmap.2(out_mat2,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main2, srtCol = srtCol, density.info = "none", cellnote=sym_mat2, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors)
      }
    } 
    
    if (interaction == TRUE & length(factor) > 1) {
      main3 = paste(main, "(Interaction)")
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat3,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main3, srtCol = srtCol, density.info = "none", cellnote=sym_mat3, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey)
      } else {
        heatmap.2(out_mat3,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main3, srtCol = srtCol, density.info = "none", cellnote=sym_mat3, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors)
      }
    }
  }
  
  if (length(factor) <= 1) {
    res = list(fac1 = out_mat1, logp1 = p_mat1)
  } else if (interaction == FALSE) {
    res = list(fac1 = out_mat1, logp1 = p_mat1, fac2 = out_mat2, logp2 = p_mat2)
  } else if (interaction == TRUE) {
    res = list(fac1 = out_mat1, logp1 = p_mat1, fac2 = out_mat2, logp2 = p_mat2, fac3 = out_mat3, logp3 = p_mat3)
  }

  return(res)
}

matrix_lm01 =  function(mat1 = matrix, 
                       mat2 = matrix, 
                       factor = "",
                       #covar = "",
                       heatmap = "FALSE", 
                       cluster = c("row", "col", "column", "both", "none"),
                       interaction = "FALSE",
                       pvalue = "FALSE", 
                       maxp = 0.05,
                       diag = "TRUE",
                       filen = "",
                       main = "",
                       pdfx = 10, pdfy = 10, rev_col = FALSE, srtCol = 90, keysize = 1, 
                       plotsymbols = TRUE, rm.other = FALSE, notecex = 1, 
                       ColSideColors = NULL , RowSideColors = NULL, minval = -10, maxval = 10,
                       flt_col_maxp = 1.0, flt_row_maxp = 1.0)
{
  
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  out_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  
  out_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  
  out_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  
  colnames(out_mat1) = colnames(mat2)
  rownames(out_mat1) = colnames(mat1)
  colnames(p_mat1) = colnames(mat2)
  rownames(p_mat1) = colnames(mat1)
  colnames(logp_mat1) = colnames(mat2)
  rownames(logp_mat1) = colnames(mat1)
  
  colnames(out_mat2) = colnames(mat2)
  rownames(out_mat2) = colnames(mat1)
  colnames(p_mat2) = colnames(mat2)
  rownames(p_mat2) = colnames(mat1)
  colnames(logp_mat2) = colnames(mat2)
  rownames(logp_mat2) = colnames(mat1)
  
  colnames(out_mat3) = colnames(mat2)
  rownames(out_mat3) = colnames(mat1)
  colnames(p_mat3) = colnames(mat2)
  rownames(p_mat3) = colnames(mat1)
  colnames(logp_mat3) = colnames(mat2)
  rownames(logp_mat3) = colnames(mat1)
  
  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      if ((diag == "TRUE") || (i != j)) {
        
        if (factor != "") {
          factest = levels(droplevels(factor[!is.na(mat2[,j])]))
        } else {
          factest = ""
        }

        if (length(factor) <= 1 | length(factest) <= 1) {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j])))
          
          out_mat1[i, j] = result[2,3]
          p_mat1[i, j] = result[2,4]
          lp1 = get_logp(result[2,4], maxp)
          logp_mat1[i, j] =lp1
          
          out_mat2[i, j] = 0
          p_mat2[i, j] = NA
          lp2 = get_logp(NA, maxp)
          logp_mat2[i, j] =lp2
          
          out_mat3[i, j] = 0
          p_mat3[i, j] = NA
          lp3 = get_logp(NA, maxp)
          logp_mat3[i, j] =lp3
          
          sym_mat1[i,j] = get_sym(lp1)
          sym_mat2[i,j] = get_sym(lp2)
          sym_mat2[i,j] = get_sym(lp3)
          
        } else if (interaction == FALSE) {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] + factor)))
          
          out_mat1[i, j] = result[2,3]
          p_mat1[i, j] = result[2,4]
          lp1 = get_logp(result[2,4], maxp)
          logp_mat1[i, j] =lp1
          
          out_mat2[i, j] = result[3,3]
          p_mat2[i, j] = result[3,4]
          lp2 = get_logp(result[3,4], maxp)
          logp_mat2[i, j] =lp2
          
          sym_mat1[i,j] = get_sym(lp1)
          sym_mat2[i,j] = get_sym(lp2)
        } else {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] * factor)))
          
          out_mat1[i, j] = result[2,3]
          p_mat1[i, j] = result[2,4]
          lp1 = get_logp(result[2,4], maxp)
          logp_mat1[i, j] =lp1
          
          out_mat2[i, j] = result[3,3]
          p_mat2[i, j] = result[3,4]
          lp2 = get_logp(result[3,4], maxp)
          logp_mat2[i, j] =lp2
          
          out_mat3[i, j] = result[4,3]
          p_mat3[i, j] = result[4,4]
          lp3 = get_logp(result[4,4], maxp)
          logp_mat3[i, j] =lp3
          
          sym_mat1[i,j] = get_sym(lp1)
          sym_mat2[i,j] = get_sym(lp2)
          sym_mat3[i,j] = get_sym(lp3)		
        }
      }
    }
  }
  
  if (pvalue == "TRUE") {
    cluster = "none"
    out_mat1 = logp_mat1
    out_mat2 = logp_mat2
    out_mat3 = logp_mat3
  }
  
  if (heatmap == "TRUE") {
    # m = out_mat1
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
      #print("both")
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
      #print("row")
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
      #print("col")
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
      #print("row")
    }
    
    if (pvalue == "TRUE") {
      if (rev_col == FALSE) {
        cols = colorpanel(100,"white","blue","red")
      } else {
        cols = colorpanel(100,"white","blue","darkblue")
      }
      xlab="-log10(P)"
    } else {
      if (rev_col == FALSE) {
        cols = colorpanel(100,"blue","white", "red")
      } else {
        cols = col=colorpanel(100,"red","white", "blue")
      }
      
      xlab="LM t-value"
    }
    
    symkey = TRUE
    
    #print(out_mat1)
    if (length(factor) > 1) {
      main1 = paste(main, "(Factor1)")
    } else {
      main1 = main
    }
    
    mx = max(out_mat1, na.rm = TRUE)
    mn = min(out_mat1, na.rm = TRUE)
    abs_mx = abs(mx)
    abs_mn = abs(mn)
    if ((mx < maxval) | (mn > minval)) {
      if (abs_mx > abs_mn) {
        mx = mx
        mn = -1 * mx
      } else {
        mx = abs_mn
        mn = mn
      }
    } else {
      mx = maxval
      mn = maxval
    }
    
    print(mx)
    print(mn)
    breaks = seq(from = mn, to = mx, by =  (mx - mn)/100)
    
    if (is.null(ColSideColors)) {
        heatmap.2(out_mat1,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                main=main1, srtCol = srtCol, density.info = "none", cellnote=sym_mat1, notecol = "white", keysize=keysize, key.xlab=xlab, 
                key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, breaks = breaks)
    } else {
        heatmap.2(out_mat1,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                main=main1, srtCol = srtCol, density.info = "none", cellnote=sym_mat1, notecol = "white", keysize=keysize, key.xlab=xlab, 
                key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors, breaks = breaks)
    }
    
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
    
    if (length(factor) > 1) {
      mx = max(out_mat2, na.rm = TRUE)
      mn = min(out_mat2, na.rm = TRUE)
      abs_mx = abs(mx)
      abs_mn = abs(mn)
      if ((mx < maxval) | (mn > minval)) {
        if (abs_mx > abs_mn) {
          mx = mx
          mn = -1 * mx
        } else {
          mx = abs_mn
          mn = mn
        }
      } else {
        mx = maxval
        mn = maxval
      }
      
      print(mx)
      print(mn)
      breaks = seq(from = mn, to = mx, by =  (mx - mn)/100)
      
      main2 = paste(main, "(Factor2)")
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat2,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main2, srtCol = srtCol, density.info = "none", cellnote=sym_mat2, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, breaks = breaks)
      } else {
        heatmap.2(out_mat2,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main2, srtCol = srtCol, density.info = "none", cellnote=sym_mat2, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors, breaks = breaks)
      }
    } 
    
    if (interaction == TRUE & length(factor) > 1) {
      mx = max(out_mat3, na.rm = TRUE)
      mn = min(out_mat3, na.rm = TRUE)
      abs_mx = abs(mx)
      abs_mn = abs(mn)
      if ((mx < maxval) | (mn > minval)) {
        if (abs_mx > abs_mn) {
          mx = mx
          mn = -1 * mx
        } else {
          mx = abs_mn
          mn = mn
        }
      } else {
        mx = maxval
        mn = maxval
      }
      
      print(mx)
      print(mn)
      breaks = seq(from = mn, to = mx, by =  (mx - mn)/100)
      
      main3 = paste(main, "(Interaction)")
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat3,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main3, srtCol = srtCol, density.info = "none", cellnote=sym_mat3, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, breaks = breaks)
      } else {
        heatmap.2(out_mat3,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main3, srtCol = srtCol, density.info = "none", cellnote=sym_mat3, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors, breaks = breaks)
      }
    }
  }
  
  if (length(factor) <= 1) {
    res = list(fac1 = out_mat1, logp1 = p_mat1)
  } else if (interaction == FALSE) {
    res = list(fac1 = out_mat1, logp1 = p_mat1, fac2 = out_mat2, logp2 = p_mat2)
  } else if (interaction == TRUE) {
    res = list(fac1 = out_mat1, logp1 = p_mat1, fac2 = out_mat2, logp2 = p_mat2, fac3 = out_mat3, logp3 = p_mat3)
  }
  
  return(res)
}

matrix_lm =  function(mat1 = matrix, 
                      mat2 = matrix, 
                      factor = "",
                      heatmap = "TRUE", 
                      cluster = c("row", "col", "column", "both", "none"),
                      interaction = "FALSE",
                      pvalue = "FALSE", 
                      maxp = 0.1,
                      diag = "TRUE",
                      filen = "",
                      main = "",
                      pdfx = 10, pdfy = 10, rev_col = FALSE, srtCol = 90, keysize = 1, 
                      plotsymbols = TRUE, rm.other = FALSE, notecex = 1, 
                      ColSideColors = NULL , RowSideColors = NULL, minval = -10, maxval = 10,
                      flt_col_maxp = 1.0, flt_row_maxp = 1.0)
{
  
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  lm_res = get_lm_matrices(mat1, mat2, factor = factor, interaction = interaction, pvalue = pvalue,
                           maxp = maxp, diag = diag)
  lm_res = shrink_lm_matrices(lm_res, flt_col_maxp = flt_col_maxp, flt_row_maxp = flt_row_maxp)
  #ncol1 = ncol(mat1)
  #ncol2 = ncol(mat2)
  #print(lm_res$fac1)
  
  out_mat1 = lm_res$fac1
  logp_mat1 = lm_res$logp1
  sym_mat1 = lm_res$sym1
  
  out_mat2 = lm_res$fac2
  logp_mat2 = lm_res$logp2
  sym_mat2 = lm_res$sym2
  
  out_mat3 = lm_res$fac3
  logp_mat3 = lm_res$logp3
  sym_mat3 = lm_res$sym3
  
  if (pvalue == "TRUE") {
    cluster = "none"
    out_mat1 = logp_mat1
    out_mat2 = logp_mat2
    out_mat3 = logp_mat3
  }
  
  if (plotsymbols == FALSE) {
    sym_mat1 = matrix(rep("", times=nrow(sym_mat1)*ncol(sym_mat1)), nrow = nrow(sym_mat1), ncol = ncol(sym_mat1))
    sym_mat2 = matrix(rep("", times=nrow(sym_mat2)*ncol(sym_mat2)), nrow = nrow(sym_mat2), ncol = ncol(sym_mat2))
    sym_mat3 = matrix(rep("", times=nrow(sym_mat3)*ncol(sym_mat3)), nrow = nrow(sym_mat3), ncol = ncol(sym_mat3))
  }
  if (heatmap == "TRUE") {
    # m = out_mat1
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
      #print("both")
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
      #print("row")
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
      #print("col")
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
      #print("row")
    }
    
    if (pvalue == "TRUE") {
      if (rev_col == FALSE) {
        cols = colorpanel(100,"white","blue","red")
      } else {
        cols = colorpanel(100,"white","blue","darkblue")
      }
      xlab="-log10(P)"
    } else {
      if (rev_col == FALSE) {
        cols = colorpanel(100,"blue","white", "red")
      } else {
        cols = col=colorpanel(100,"red","white", "blue")
      }
      
      xlab="LM t-value"
    }
    
    symkey = TRUE
    
    if (length(factor) > 1) {
      main1 = paste(main, "(Factor1)")
    } else {
      main1 = main
    }
    
    mx = max(out_mat1, na.rm = TRUE)
    mn = min(out_mat1, na.rm = TRUE)
 #   print(mn)
 #   print(mx)
    abs_mx = abs(mx)
    abs_mn = abs(mn)
    if ((mn < minval) | (mx > maxval)) {
      mx = maxval
      mn = minval
      
    } else {
      if (abs_mx > abs_mn) {
        mx = mx
        mn = -1 * mx
      } else {
        mx = abs_mn
        mn = mn
      }
    }
    breaks = seq(from = mn, to = mx, by =  (mx - mn)/100)
#    print(mn)
#    print(mx)
    cat("\n")
    if (is.null(ColSideColors)) {
      heatmap.2(out_mat1,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                main=main1, srtCol = srtCol, density.info = "none", cellnote=sym_mat1, notecol = "white", keysize=keysize, key.xlab=xlab, 
                key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, breaks = breaks)
    } else {
      heatmap.2(out_mat1,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                main=main1, srtCol = srtCol, density.info = "none", cellnote=sym_mat1, notecol = "white", keysize=keysize, key.xlab=xlab, 
                key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors, breaks = breaks)
    }
    
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
    
    if (length(factor) > 1) {
      mx = max(out_mat2, na.rm = TRUE)
      mn = min(out_mat2, na.rm = TRUE)
 #     print(mx)
 #     print(mn)
      
      abs_mx = abs(mx)
      abs_mn = abs(mn)
      if ((mn < minval) | (mx > maxval)) {
        mx = maxval
        mn = minval
      } else {
        if (abs_mx > abs_mn) {
          mx = mx
          mn = -1 * mx
        } else {
          mx = abs_mn
          mn = mn
        }
      }
      
#      print(mx)
#      print(mn)
      breaks = seq(from = mn, to = mx, by =  (mx - mn)/100)
      
      #par(mar=c(6.0, 6.0, 4.0, 1.0))
      main2 = paste(main, "(Factor2)")
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat2,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main2, srtCol = srtCol, density.info = "none", cellnote=sym_mat2, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, breaks = breaks)
      } else {
        heatmap.2(out_mat2,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main2, srtCol = srtCol, density.info = "none", cellnote=sym_mat2, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors, breaks = breaks)
      }
    } 
    cat("\n")
    if (interaction == TRUE & length(factor) > 1) {
      mx = max(out_mat3, na.rm = TRUE)
      mn = min(out_mat3, na.rm = TRUE)
      #print(mx)
      #print(mn)
      
      abs_mx = abs(mx)
      abs_mn = abs(mn)
      if ((mn < minval) | (mx > maxval)) {
        mx = maxval
        mn = minval
       
      } else {
        if (abs_mx > abs_mn) {
          mx = mx
          mn = -1 * mx
        } else {
          mx = abs_mn
          mn = mn
        }
      }
      
      #print(mx)
      #print(mn)
      breaks = seq(from = mn, to = mx, by =  (mx - mn)/100)
      
      main3 = paste(main, "(Interaction)")
      if (is.null(ColSideColors)) {
        heatmap.2(out_mat3,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main3, srtCol = srtCol, density.info = "none", cellnote=sym_mat3, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, breaks = breaks)
      } else {
        heatmap.2(out_mat3,col=cols,trace="none", Colv=colv, Rowv=rowv, 
                  main=main3, srtCol = srtCol, density.info = "none", cellnote=sym_mat3, notecol = "white", keysize=keysize, key.xlab=xlab, 
                  key.title = NA, notecex = notecex, symkey = symkey, symbreaks = symkey, ColSideColors = ColSideColors, breaks = breaks)
      }
    }
  }
  
  if (length(factor) <= 1) {
    res = list(fac1 = out_mat1, logp1 = logp_mat1)
  } else if (interaction == FALSE) {
    res = list(fac1 = out_mat1, logp1 = logp_mat1, fac2 = out_mat2, logp2 = logp_mat2)
  } else if (interaction == TRUE) {
    res = list(fac1 = out_mat1, logp1 = logp_mat1, fac2 = out_mat2, logp2 = logp_mat2, fac3 = out_mat3, logp3 = logp_mat3)
  }
  
  return(res)
}

shrink_lm_matrices =  function(lm_res, flt_col_maxp = 1.0, flt_row_maxp = 1.0)
{
  #res = list(fac1 = out_mat1, logp1 = logp_mat1, sym1 = sym_mat1, fac2 = out_mat2, logp2 = logp_mat2, sym2 = sym_mat2, fac3 = out_mat3, logp3 = logp_mat3, sym3 = sym_mat3)
  
  out_mat1 = lm_res$fac1
  #print(out_mat1)
  cat("Initial NCOL: ", ncol(out_mat1), "\n")
  cat("Initial NROW: ", nrow(out_mat1), "\n")
  
  logp_mat1 = lm_res$logp1
  sym_mat1 = lm_res$sym1
  
  out_mat2 = lm_res$fac2
  logp_mat2 = lm_res$logp2
  sym_mat2 = lm_res$sym2
  
  out_mat3 = lm_res$fac3
  logp_mat3 = lm_res$logp3
  sym_mat3 = lm_res$sym3
  
  rows = nrow(out_mat1)
  inc_rows = rep(TRUE, rows)
  cols = ncol(out_mat1)
  inc_cols = rep(TRUE, cols)
  
  log_col_maxp = -1 * log10(flt_col_maxp)
  log_row_maxp = -1 * log10(flt_row_maxp)
#  print(log_col_maxp)
#  print(log_row_maxp)
  
  for (i in 1:rows) {
    max1 = max(logp_mat1[i, ], na.rm = TRUE)
    max2 = max(logp_mat2[i, ], na.rm = TRUE)
    max3 = max(logp_mat2[i, ], na.rm = TRUE)
#    cat(max1, max2, max3, "\n")
    if (max(cbind(max1,max2,max3)) < log_row_maxp) inc_rows[i] = FALSE
    
    if ((sum(!is.na(out_mat1[i, ])) < 1) | (sum(!is.na(out_mat2[i, ])) < 1) | (sum(!is.na(out_mat3[i, ])) < 1)) inc_rows[i] = FALSE
  }
#  print( sum(inc_rows))
  
  #Only alter matrices if at least two rows will be retained
  if (sum(inc_rows) > 1) {
    out_mat1 = out_mat1[inc_rows == TRUE, ]
    logp_mat1 = logp_mat1[inc_rows == TRUE, ]
    sym_mat1 = sym_mat1[inc_rows == TRUE, ]
    
    out_mat2 = out_mat2[inc_rows == TRUE, ]
    logp_mat2 = logp_mat2[inc_rows == TRUE, ]
    sym_mat2 = sym_mat2[inc_rows == TRUE, ]
    
    out_mat3 = out_mat3[inc_rows == TRUE, ]
    logp_mat3 = logp_mat3[inc_rows == TRUE, ]
    sym_mat3 = sym_mat3[inc_rows == TRUE, ]
    
  }
  
  for (i in 1:cols) {
    max1 = max(logp_mat1[, i], na.rm = TRUE)
    max2 = max(logp_mat2[, i], na.rm = TRUE)
    max3 = max(logp_mat2[, i], na.rm = TRUE)
    if (max(cbind(max1,max2,max3)) < log_col_maxp) inc_cols[i] = FALSE
    
    if ((sum(!is.na(out_mat1[, i])) < 1) | (sum(!is.na(out_mat2[, i])) < 1) | (sum(!is.na(out_mat3[, i])) < 1)) inc_rows[i] = FALSE
    
#    cat(max1, " ", max2, " ", max3, " ", max(cbind(max1,max2,max3)), "\n")
    
  }
  #Only alter matrices if at least two columns will be retained
  if (sum(inc_cols) > 1) {
    out_mat1 = out_mat1[, inc_cols == TRUE]
    logp_mat1 = logp_mat1[, inc_cols == TRUE]
    sym_mat1 = sym_mat1[, inc_cols == TRUE]
    
    out_mat2 = out_mat2[, inc_cols == TRUE]
    logp_mat2 = logp_mat2[, inc_cols == TRUE]
    sym_mat2 = sym_mat2[, inc_cols == TRUE]
    
    out_mat3 = out_mat3[, inc_cols == TRUE]
    logp_mat3 = logp_mat3[, inc_cols == TRUE]
    sym_mat3 = sym_mat3[, inc_cols == TRUE]
  }
  
  cat("Final NCOL: ", ncol(out_mat1), "\n")
  cat("Final NROW: ", nrow(out_mat1), "\n")
  
  res = list(fac1 = out_mat1, logp1 = logp_mat1, sym1 = sym_mat1, fac2 = out_mat2, logp2 = logp_mat2, sym2 = sym_mat2, fac3 = out_mat3, logp3 = logp_mat3, sym3 = sym_mat3)
  print("End shrink_lm_matrices")
  return(res)
}

get_lm_matrices =  function(mat1 = matrix, 
                            mat2 = matrix, 
                            factor = "",
                            interaction = "FALSE",
                            pvalue = "FALSE", 
                            maxp = 0.05,
                            diag = "TRUE",
                            flt_col_maxp = 1.0, flt_row_maxp = 1.0)
{
  
  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  out_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat1 = matrix(nrow = ncol1, ncol = ncol2)
  
  out_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat2 = matrix(nrow = ncol1, ncol = ncol2)
  
  out_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  p_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  logp_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat3 = matrix(nrow = ncol1, ncol = ncol2)
  
  colnames(out_mat1) = colnames(mat2)
  rownames(out_mat1) = colnames(mat1)
  colnames(p_mat1) = colnames(mat2)
  rownames(p_mat1) = colnames(mat1)
  colnames(logp_mat1) = colnames(mat2)
  rownames(logp_mat1) = colnames(mat1)
  
  colnames(out_mat2) = colnames(mat2)
  rownames(out_mat2) = colnames(mat1)
  colnames(p_mat2) = colnames(mat2)
  rownames(p_mat2) = colnames(mat1)
  colnames(logp_mat2) = colnames(mat2)
  rownames(logp_mat2) = colnames(mat1)
  
  colnames(out_mat3) = colnames(mat2)
  rownames(out_mat3) = colnames(mat1)
  colnames(p_mat3) = colnames(mat2)
  rownames(p_mat3) = colnames(mat1)
  colnames(logp_mat3) = colnames(mat2)
  rownames(logp_mat3) = colnames(mat1)
  
  for (i in 1:ncol1) {
   # if (sum(!is.na(mat1[,i])) < 1) next
      for (j in 1:ncol2) {
        #cat(i,j,"\n")
        #print(mat1[,i])
        #print(mat2[,j])
        #print(sum(!is.na(mat1[,i])))
        #print(sum(!is.na(mat2[,j])))
        nas = !is.na(mat1[,i]) & !is.na(mat2[,j])
        #print(sum(nas))
        #print(nas)
        if (((sum(!is.na(mat2[,j])) < 1) | (sum(!is.na(mat1[,i])) < 1)) | (sum(nas) <= 2)) {
          out_mat1[i, j] = NA
          out_mat2[i, j] = NA
          out_mat3[i, j] = NA
          logp_mat1[i, j] = get_logp(NA, maxp)
          logp_mat2[i, j] = get_logp(NA, maxp)
          logp_mat3[i, j] = get_logp(NA, maxp)
 #         cat("Skipping: ", i,j, "\n")
          next
        }
        
        if ((diag == "TRUE") | (i != j)) {
          if (factor != "") {
            factest = levels(droplevels(factor[!is.na(mat2[,j])]))
          } else {
            factest = ""
          }
        
          if (length(factor) <= 1 | length(factest) <= 1) {
            result = coefficients(summary(lm(mat1[,i] ~ mat2[,j])))
          
            out_mat1[i, j] = result[2,3]
            p_mat1[i, j] = result[2,4]
            lp1 = get_logp(result[2,4], maxp)
            logp_mat1[i, j] =lp1
          
            out_mat2[i, j] = 0
            p_mat2[i, j] = NA
            lp2 = get_logp(NA, maxp)
            logp_mat2[i, j] =lp2
          
            out_mat3[i, j] = 0
            p_mat3[i, j] = NA
            lp3 = get_logp(NA, maxp)
            logp_mat3[i, j] =lp3
          
            sym_mat1[i,j] = get_sym(lp1)
            sym_mat2[i,j] = get_sym(lp2)
            sym_mat2[i,j] = get_sym(lp3)
          
          } else if (interaction == FALSE) {
            
            result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] + factor)))
          
            out_mat1[i, j] = result[2,3]
            p_mat1[i, j] = result[2,4]
            lp1 = get_logp(result[2,4], maxp)
            logp_mat1[i, j] =lp1
          
            out_mat2[i, j] = result[3,3]
            p_mat2[i, j] = result[3,4]
            lp2 = get_logp(result[3,4], maxp)
            logp_mat2[i, j] =lp2
          
            sym_mat1[i,j] = get_sym(lp1)
            sym_mat2[i,j] = get_sym(lp2)
          } else {
           # print(mat1[,i])
            #print(mat2[,j])

            result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] * factor)))
          
            out_mat1[i, j] = result[2,3]
            p_mat1[i, j] = result[2,4]
            lp1 = get_logp(result[2,4], maxp)
            logp_mat1[i, j] =lp1
          
            out_mat2[i, j] = result[3,3]
            p_mat2[i, j] = result[3,4]
            lp2 = get_logp(result[3,4], maxp)
            logp_mat2[i, j] =lp2
          
            out_mat3[i, j] = result[4,3]
            p_mat3[i, j] = result[4,4]
            lp3 = get_logp(result[4,4], maxp)
            logp_mat3[i, j] =lp3
          
            sym_mat1[i,j] = get_sym(lp1)
            sym_mat2[i,j] = get_sym(lp2)
            sym_mat3[i,j] = get_sym(lp3)		
          }
        }
    }
  }
  res = list(fac1 = out_mat1, logp1 = logp_mat1, sym1 = sym_mat1, fac2 = out_mat2, logp2 = logp_mat2, sym2 = sym_mat2, fac3 = out_mat3, logp3 = logp_mat3, sym3 = sym_mat3)
  #print(out_mat1)
  print("End get_lm_matrices")
  return(res)
}


get_logp =  function(p = p, maxp = maxp)
{
  if (is.na(p) == TRUE) {
    p = 1.0
  } else if (p < 0.0000001) {
    p = 0.000001
  } 
  lp = -1*log10(p)
  if (lp < -1*log10(maxp)) lp = 0
  return(lp)
}

get_sym = function(lp = lp)
{
  sym = ''
  if (lp >= -log10(0.001)) {
    sym = '#'
  } else if (lp >= -log10(0.01)) {
    sym = '+'
  } else if (lp >= -log10(0.05)) {
    sym = '*'
  } else if (lp >= -log10(0.1)) {
    sym = 'o'
  }
  return(sym)
}

old_matrix_lm =  function(mat1 = matrix, 
                        mat2 = matrix, 
                        factor = "",
                        #covar = "",
                        heatmap = "FALSE", 
                        cluster = c("row", "col", "column", "both", "none"),
                        pvalue = "FALSE", 
                        maxp = 0.05,
                        diag = "TRUE",
                        filen = "",
                        main = "",
                        pdfx = 10, pdfy = 10, rev_col = FALSE, srtCol = 90, keysize = 1, 
                        plotsymbols = TRUE, rm.other = FALSE, notecex = 1, ColSideColors = "" )
{

  if (nrow(mat1) != nrow(mat2)) stop("Both matrices must have the same number of rows")
  
  if (rm.other == TRUE) {
    mat1 = mat1[, colnames(mat1) != "Other"]
    mat2 = mat2[, colnames(mat2) != "Other"]
  }
  
  ncol1 = ncol(mat1)
  ncol2 = ncol(mat2)
  out_mat = matrix(nrow = ncol1, ncol = ncol2)
  p_mat = matrix(nrow = ncol1, ncol = ncol2)
  sym_mat = matrix(nrow = ncol1, ncol = ncol2)
  colnames(out_mat) = colnames(mat2)
  rownames(out_mat) = colnames(mat1)
  colnames(p_mat) = colnames(mat2)
  rownames(p_mat) = colnames(mat1)

  for (i in 1:ncol1) {
    for (j in 1:ncol2) {
      if ((diag == "TRUE") || (i != j)) {
        
        if (length(factor) <= 1) {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j])))
        } else {
          result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] + factor)))
          #result = coefficients(summary(lm(mat1[,i] ~ mat2[,j] + factor + covar)))
        }
        
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
        
        out_mat[i, j] = result[2,3]

        p = result[2,4]
        if (is.na(p) == TRUE) {
          p = 1.0
        } else if (p < 0.0000001) {
            p = 0.000001
        } 
        lp = -1*log10(p)
        #lp = p
        if (lp < -1*log10(maxp)) lp = 0
        p_mat[i, j] = lp
    
        if (lp >= -log10(0.01) & plotsymbols == TRUE) {
          sym_mat[i, j] = '#'
        } else if (lp >= -log10(0.05) & plotsymbols == TRUE) {
          sym_mat[i, j] = '+'
        }
      }
    }
  }
  
  if (pvalue == "TRUE") {
      cluster = "none"
      out_mat = p_mat
  }
  
  if (heatmap == "TRUE")
  {
    m = out_mat
    colv = FALSE
    rowv = FALSE
    if (cluster == "both") {
      colv = TRUE
      rowv = TRUE
      #print("both")
    } else if (cluster == "row") {
      colv = FALSE
      rowv = TRUE
      #print("row")
    } else if ((cluster == "column") || (cluster == "col")) {
      colv = TRUE
      rowv = FALSE
      #print("col")
    } else if (cluster == "none") {
      colv = FALSE
      rowv = FALSE
      #print("row")
    }
    
    if (ColSideColors == "") ColSideColors = rep("white",ncol(m))
    
    if (pvalue == "TRUE") {
       if (rev_col == FALSE) {
        heatmap.2(m,col=colorpanel(100,"white","blue","red"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE, ColSideColors = ColSideColors)
      } else {
        heatmap.2(m,col=colorpanel(100,"white","blue","darkblue"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-log10(P)", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE, ColSideColors = ColSideColors)
      }
    } else {
      if (rev_col == FALSE) {
        heatmap.2(m,col=colorpanel(100,"blue","white", "red"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="LM t-value", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE, ColSideColors = ColSideColors)
      } else {
        heatmap.2(m,col=colorpanel(100,"red","white", "blue"),trace="none", Colv=colv, Rowv=rowv, 
                  main=main, srtCol = srtCol, density.info = "none", cellnote=sym_mat, notecol = "white", keysize=keysize, key.xlab="-LM t-value", 
                  key.title = NA, notecex = notecex, symkey = FALSE, symbreaks = FALSE, ColSideColors = ColSideColors)
      }
    }
    
    #heatmap.2(m,trace="none", Colv=colv, Rowv=rowv)
    if (filen != "") {
      dev.copy(pdf, filen, width = pdfx, height = pdfy)
      dev.off()
    }
  }
  
  return(out_mat)
}