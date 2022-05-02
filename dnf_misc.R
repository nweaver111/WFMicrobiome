#************************
#
#
# v.0.1.0  
# daniel.frank@ucdenver.edu
# 16 Sept 2018
#
#
#************************

library(Hmisc)  #For drawing minor tick marks


read.otutable <- function(fn, order = TRUE, rm.prev = TRUE, otu.nms="OTU_Name", transform = TRUE, meta = "", drop.zeros = TRUE, explicet = FALSE, sep = '\t') {
  cat("Reading file:", fn, "\n")
  f_cts <- read.table(fn, header = TRUE, row.names=otu.nms, na.strings = c("na","NA","nd","n.a.","n.d."), check.names = FALSE, sep = sep)

  if (rm.prev == TRUE && explicet == FALSE) f_cts <- f_cts[,-1]  #Remove prevalence
  if (explicet == TRUE) f_cts <- f_cts[-1, ]  #Remove root row
  f_cts <- t(f_cts)  # transpose
  
  orignms = colnames(f_cts)
  new_otu_nms = otu.shortnames(orignms)
  colnames(f_cts) <- new_otu_nms$short
  
  if (order == TRUE) f_cts = f_cts[order(rownames(f_cts), decreasing = FALSE), ]
  f_trans = list(cts=f_cts, per="", clr="", bin="", nms="")
  if (transform == TRUE) f_trans = mul_transform(f_cts)
  obj = (list(cts = f_trans$cts, per = f_trans$per, clr = f_trans$clr, bin = f_trans$bin, nms = f_trans$nms, meta = meta, phynms = new_otu_nms$phynms, orignms = orignms))
  if (drop.zeros == TRUE) obj = otu.drop_zero_counts(obj)
  
  if (length(rownames(meta)) != length(rownames(f_cts))) stop("Metadata and OTU tables differ in number of subjects!")

  mis_count = 0
  for (i in 1:length(rownames(meta))) {
    
    if (rownames(meta)[i] != rownames(f_cts)[i]) {
        mis_count = mis_count + 1
     }
  }
  if (mis_count > 0) {
    stop("Subject names in metadata and OTU tables do not match!")
  } else {
    cat("Found", length(rownames(f_cts)), "subjects in OTU file.\n" )
    cat("Found", length(colnames(f_cts)), "OTUs in OTU file.\n" )
  }
  
  if (anyDuplicated(rownames(f_cts))) stop("OTU file contains duplicate subject/library names!")
      
  return(obj)
}

read.metatable <- function(fn, order = TRUE, lib.nms="Lib") {
  cat("Reading file:", fn, "\n")
  f_meta <- read.delim(fn, header = TRUE, na.strings = c("na","NA","nd","n.a.","n.d.", "N/A", "n/a"), sep = "\t")
  #print(f_meta[,lib.nms])
  rownames(f_meta) <- f_meta[,lib.nms]
  if (order == TRUE) f_meta = f_meta[order(rownames(f_meta), decreasing = FALSE), ]
  
  cat("Found", length(rownames(f_meta)), "subjects in file.\n" )
  cat("Found", length(colnames(f_meta)), "variables in file.\n" )
  if (anyDuplicated(rownames(f_meta))) stop("Metadata file contains duplicate subject/library names!")
  
  return(f_meta)
}

otu.summary <- function(obj) {
  cat("\n")
  cat("OTU Table: ",deparse(substitute(obj)),"\n")
  cat(" Subjects: ", nrow(obj$cts), "\n")
  cat("     OTUs: ", ncol(obj$cts), "\n")
  cat("    Reads: ", sum(obj$cts), "\n")
  cat("Meta Vars: ", ncol(obj$meta), "\n\n")
  
}

otu.dim <- function(obj) {
  l = list(nsubj = nrow(obj$cts), notu = ncol(obj$cts))
  return(l)
}

otu.nsubj <- function(obj) {
  return(nrow(obj$cts))
}

otu.notu <- function(obj) {
  return(ncol(obj$cts))
}

otu.subj_names <- function(obj) {
  return(rownames(obj$cts))
}

otu.summarize_read_cts <- function(obj, mask = "", factor ="", logtrans = FALSE, boxplot = FALSE) {
  RND = 1
  
  if ((mask != "" && factor != "") && (length(mask) != length(factor))) stop("Mask and factor must have same lengths.\n")
  cts = otu.read_cts(obj)
  if (logtrans == TRUE) cts = log10(cts)
  if (length(mask) <= 1) mask = !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  #if (length(factor) <= 1) factor = !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  if (length(factor) <= 1) factor = rep("ALL", otu.nsubj(obj))
  
  Reads = otu.read_cts(obj)
  cts = cts[mask]
  factor = factor[mask]
  
  n  = aggregate(cts, by = list(factor), FUN = length)
  min = aggregate(cts, by = list(factor), FUN = min)
  max = aggregate(cts, by = list(factor), FUN = max)
  mn = aggregate(cts, by = list(factor), FUN = mean)
  sd = aggregate(cts, by = list(factor), FUN = sd)
  md = aggregate(cts, by = list(factor), FUN = median)
  q1 = aggregate(cts, by = list(factor), FUN = quantile, probs=0.25)
  q3 = aggregate(cts, by = list(factor), FUN = quantile, probs=0.75)
  
  if (logtrans == TRUE) rnd = RND
  else rnd = 0

  mat = cbind(n$x, min$x, max$x, mn$x, sd$x, md$x, q1$x, q3$x)
  mat = round(mat, rnd)
  
  res = data.frame(mat, row.names = min$Group.1)
  colnames(res) = c("N", "min", "max", "mean", "sd", "med", "q1", "q3")
  
  st1 = aggregate(cts, by = list(factor), median_iqr_str, rnd = rnd)
  res["Med (IQR)"] = st1$x
  st2 = aggregate(cts, by = list(factor), mean_sd_str, rnd = rnd)
  res["Mean (SD)"] = st2$x
  
  if (boxplot == TRUE) {
    if (logtrans == TRUE) ylab = "log10(Reads)"
    else ylab = "Reads"
    opar = par()
    par(mfrow=c(1,1))
    boxplot(cts ~ factor, main="Read Counts", ylab=ylab)
    par(opar)
  }
  return(res)
}

otu.read_cts <- function(obj) {
  cts = apply(obj$cts, MARGIN = 1, sum)
  return(cts)
}

otu.subset_libs <- function(obj, mask = mask, drop.zeros = TRUE, vrb = TRUE) {
  mask = mask & !is.na(mask) #Account for NAs
  
  cts = (obj$cts)
  cts = cts[mask,]
  meta = (obj$meta)
  meta = meta[mask,]
  meta = droplevels(meta)
  obj$cts = cts
  obj$meta = meta
  
  if (drop.zeros == TRUE) obj = otu.drop_zero_counts(obj, vrb)
  obj = otu.mul_transform(obj)
  return(obj)
}

otu.group_ra <- function(obj, mask = "", factor ="", ra_cutoff = 0.0) {
  if ((mask != "" && factor != "") && (length(mask) != length(factor))) stop("Mask and factor must have same lengths.\n")
  if (length(mask) <= 1) mask = !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  
  #print(length(factor))
  #if (length(factor) <= 1) factor = rep("ALL", otu.nsubj(obj))
  
  if (length(factor) <= 1) {
    if (tolower(factor) == "all") factor = rep("ALL", otu.nsubj(obj))
    else factor = otu.subj_names(obj)
  }

  #print(factor)
  
  if (ra_cutoff < 0.0 || ra_cutoff > 100.0) stop("RA cutoff must be 0.0 <= cutoff <= 100.0.\n")
  per = (obj$per)[mask,]
  factor = factor[mask]
  
  
  mn = aggregate(100*per, by = list(factor), FUN = mean)
  md = aggregate(100*per, by = list(factor), FUN = median)

  rownames(mn) = mn[,1]
  rownames(md) = md[,1]
  
  mn = mn[,-1]
  md = md[,-1]
  other = numeric(nrow(mn))
 

  if (ra_cutoff > 0.0) {
    mx = apply(mn, MARGIN = 2, max)
    
    orig = mn
    mn = mn[,(mx >= ra_cutoff)]
    md = md[,(mx >= ra_cutoff)]
    
    x = as.matrix(orig[,(mx < ra_cutoff)])
    
    if (ncol(x) == 1) {
      rownames(x) = rownames(orig)
      colnames(x) = (colnames(orig))[(mx < ra_cutoff)]
    }
    other = apply(x, MARGIN = 1, sum)
  }
  return(list(mean=mn, median=md, other=other))
}

otu.summarize_ra <- function(obj, mask = "") {
  if (length(mask) <= 1) mask = !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  per = (obj$per)[mask,]
  notus = ncol(per)

  factor = rep("ALL", otu.nsubj(obj))
  
  mns = vector(mode = "numeric", length = notus)
  mds = vector(mode = "numeric", length = notus)
  mins = vector(mode = "numeric", length = notus)
  maxs = vector(mode = "numeric", length = notus)
  sds = vector(mode = "numeric", length = notus)
  q1s = vector(mode = "numeric", length = notus)
  q3s = vector(mode = "numeric", length = notus)
  st1s = vector(mode = "numeric", length = notus)
  st2s = vector(mode = "numeric", length = notus)
  rnd = 4
  
  for (i in 1:notus) {
    mns[i] = (aggregate(100*per[,i], by = list(factor), FUN = mean))$x
    mds[i] = (aggregate(100*per[,i], by = list(factor), FUN = median))$x
    mins[i] = (aggregate(100*per[,i], by = list(factor), FUN = min))$x
    maxs[i] = (aggregate(100*per[,i], by = list(factor), FUN = max))$x
    sds[i] = (aggregate(100*per[,i], by = list(factor), FUN = sd))$x
    q1s[i] = (aggregate(100*per[,i], by = list(factor), FUN = quantile, probs=0.25))$x
    q3s[i] = (aggregate(100*per[,i], by = list(factor), FUN = quantile, probs=0.75))$x
    st1s[i] = (aggregate(100*per[,i], by = list(factor),  median_iqr_str, rnd = rnd))$x
    st2s[i] = (aggregate(100*per[,i], by = list(factor), mean_sd_str, rnd = rnd))$x
  }
  
  mns = round(mns,rnd)
  mds = round(mds,rnd)
  mins = round(mins,rnd)
  maxs = round(maxs,rnd)
  sds = round(sds,rnd)
  q1s = round(q1s,rnd)
  q3s = round(q3s,rnd)

  mat = cbind(mins,maxs,mns,sds,mds,q1s,q3s,st1s,st2s)
  res = data.frame(mat,row.names=colnames(per))
  
  colnames(res) = c("min", "max", "mean", "sd", "med", "q1", "q3","Med (IQR)","Mean (SD)")
  return(res)
}


otu.plot_ra_bar <- function(obj, mask = "", factor ="", ra_cutoff = 0.0, 
                            max_cutoff = 100, horiz = FALSE, pal = "rainbow", 
                            revcol = FALSE, mix = FALSE, results = FALSE, las = 1, 
                            rm.other = FALSE, main="") {
  
  grp = otu.group_ra(obj, mask, factor, ra_cutoff)
  col_count = ncol(grp$mean)

  if ("Other" %in% colnames(grp$mean)) col_count = col_count -1
  
  if (pal == "rich") cols = rich.colors(col_count)
  else if (pal == "rainbow") cols = rainbow(col_count)
  else if (pal == "terrain") cols = terrain.colors(col_count)
  else if (pal == "topo") cols = topo.colors(col_count)
  
  if (revcol == TRUE) cols = rev(cols)
  
  if (mix == TRUE) {
    mid = round(length(cols)/2)
    a = cols[1:mid]
    b = cols[(mid+1):length(cols)]
    b = rev(b)
    mlab <- min(length(a), length(b)) 
    seqmlab <- seq(length=mlab)
    newcols <- c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab]) 
    cols = newcols
  }
  

  if ("Other" %in% colnames(grp$mean)) {
    o = (grp$mean)[,"Other"]
    other_mask = colnames(grp$mean) == "Other"
    origrows = rownames(grp$mean)
    origcols = colnames(grp$mean)
    grp$mean = as.matrix((grp$mean)[,!other_mask])
    
    #If only one OTU is selected, results are returned as a vector
    #Must coerce into matrix and add row and col names
    if (ncol(grp$mean) == 1) {
      newcols = origcols[!other_mask]
      #print(newcols)
      colnames(grp$mean) = newcols
      rownames(grp$mean) = origrows
    }
    o = o + grp$other
    o = as.matrix(o)
    colnames(o) = "Other"
    if (rm.other == FALSE) {
      grp$mean = cbind(grp$mean, o)
    }
  } else {
    o = as.matrix(grp$other)
    colnames(o) = "Other"
    if (rm.other == FALSE) {
      grp$mean = cbind(grp$mean, o)
    }
  }
  
  #Note: reversing legend names and colors synchs them with order of colors in barplot
  legnames = rev(colnames(grp$mean))
  ptcols = rev(cols)
  
  if (rm.other == FALSE) {
    cols = c(cols, "#FFFFFF")
    ptcols = c("#FFFFFF", ptcols)
  }
  
  m <- matrix(c(1, 1,1,1,1,1,1, 2,2,2,2,3), nrow = 1, ncol = 12, byrow=TRUE)
  n = layout(m)
  #layout.show(n)
  
  if (horiz == TRUE) {
    xl="rRNA Abundance (%)"
    yl=""
  } else {
    yl="rRNA Abundance (%)"
    xl=""
  }
  
  opar = par()
  par(mar=c(6.0, 6.0, 4.0, 1.0))
  if (main == "") main = deparse(substitute(factor))
  
  mx = 0.0
  grps = nrow(grp$mean)
  
  if (!is.null(grps)) {
    for (g in 1:grps) {
      if (sum((grp$mean)[g, ]) > mx) mx = sum((grp$mean)[g, ])
    }
  } else {
    mx = sum(grp$mean)
  }
  
  barplot(t(grp$mean), horiz = horiz, col = cols, 
          las=las, xlab=xl, ylab=yl, 
          cex.names=1.5, cex.lab = 1.5, 
          main=main, cex.main=2,ylim=range(pretty(c(0, mx))))
  par(mar=c(0, 1.5, 0, 0))
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
  legend("left", legend = legnames, horiz = FALSE, xpd=NA, cex= 1, pch=22, pt.bg = ptcols)
  par(opar)
  if (results == TRUE) return(grp)
}

otu.plot_ra_heatmap <- function(obj, mask = "", factor ="", ra_cutoff = 0.0, horiz = FALSE, pal = "rich", results = FALSE, main="") {
  grp = otu.group_ra(obj, mask, factor, ra_cutoff)
  
  col_count = ncol(grp$mean)
  if ("Other" %in% colnames(grp$mean)) col_count = col_count -1
  
  if (pal == "rich") cols = rich.colors(col_count)
  else if (pal == "rainbow") cols = rainbow(col_count)
  
  if ("Other" %in% colnames(grp$mean)) {
    o = (grp$mean)[,"Other"]
    other_mask = colnames(grp$mean) == "Other"
    grp$mean = (grp$mean)[,!other_mask]
    o = o + grp$other
    o = as.matrix(o)
    colnames(o) = "Other"
    grp$mean = cbind(grp$mean, o)
  } else {
    o = as.matrix(grp$other)
    colnames(o) = "Other"
    grp$mean = cbind(grp$mean, o)
  }
  if (horiz == TRUE) mat = t(grp$mean)
  else mat = as.matrix(grp$mean)
  #heatmap.2(mat,col=colorpanel(100,"white","blue","red"),trace="none", dendrogram="none", Rowv=FALSE, Colv = FALSE, main = main,lmat=rbind( c(0, 3, 4), c(2,1,0) ), lwid=c(1.5, 4, 2 ), srtCol=45)
  if (main == "") main = deparse(substitute(factor))
  opar=par()
  par(oma=c(0,0,2,8))
  heatmap.2(mat,col=colorpanel(100,"white","blue","red"),trace="none", dendrogram="none", Rowv=FALSE, Colv = FALSE, main = main, srtCol=45)
  par(opar)
  if (results == TRUE) return(grp)
}

otu.filter_otus <- function(obj, prevcut = 0.0, racut = 0.0, unc = FALSE, bact = FALSE, transform = TRUE, vrb = TRUE) {
  if (prevcut > 0.0) obj = otu.filter_prev(obj, prevcut, transform=FALSE, vrb = vrb)
  if (racut > 0.0) obj = otu.filter_ra(obj, racut, transform=FALSE, vrb = vrb)
  if (unc == TRUE) obj = otu.filter_unc(obj, transform = FALSE, vrb = vrb)
  if (bact == TRUE) obj = otu.filter_bact(obj, transform = FALSE, vrb = vrb)
  if (transform == TRUE) obj = otu.mul_transform(obj)
  return(obj)
}

otu.filter_prev <- function(obj, cutoff = cutoff, transform = TRUE, vrb = TRUE) {
  if ((cutoff > 100) | (cutoff < 0)) stop("Cutoff must be between 0 and 100 (%)")
  if (cutoff == 0) {
    if (vrb == TRUE) cat("Object unchanged since prevalence cutoff set to 0%")
    return(obj)
  }
  starting_cts = sum(obj$cts, na.rm = TRUE)
  if (vrb == TRUE) cat("Prevalence cutoff:",cutoff,"\b% (i.e., at least",cutoff,"\b% of libaries must be represented to keep OTU)\n")
  cts = obj$cts
  prevs = apply(cts, MARGIN = 2, otu.prev)
  below_cutoff = sum(prevs < cutoff, na.rm = TRUE)
  if (below_cutoff == ncol(obj$cts)) stop("All OTUs would be excluded with this cutoff. Halting operationA")
  
  #if (sum(exclude) > 0) {
  if (below_cutoff > 0) {
    include = cts[,prevs >= cutoff]
    exclude = cts[,prevs < cutoff]
    
    obj$orignms = obj$orignms[prevs >= cutoff]
    obj$phynms = obj$phynms[prevs >= cutoff]
    
    if (class(exclude) == "matrix") other = apply(exclude, MARGIN = 1, sum)
    else if (class(exclude) == "integer") other = exclude
    
    if ("Other" %in% colnames(include)) {
      include[ ,"Other"] <- include[ ,"Other"] + other
      cat("Found \"Other\" category in input data.\n")
    } else {
      other = as.matrix(other)
      colnames(other) <- "Other"

      obj$orignms <- append(obj$orignms, "Other")
      obj$phynms <- append(obj$phynms, "Other")
      include <- cbind(include, other)
      cat("Created new \"Other\" category.\n")
    }
    
    if (sum(include, na.rm = TRUE) != sum(obj$cts, na.rm = TRUE)) stop("Error in filtering OTU table!")
  
    per = clr = bin = ""
    if (transform == TRUE) {
      per = per_transform(include)
      clr = clr_transform(include)
      bin = bin_transform(include)
    }
    if (vrb == TRUE) {
      cat("Collapsed",below_cutoff,"OTUs to \"Other\" otu category.\n")
      cat("Converted",sum(exclude),"counts to \"Other\" otu category.\n")
      cat("Remaining OTUS:",ncol(include)," (Including \"Other\").\n\n")
    }
    obj$cts = include
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(include)
    
    ending_cts = sum(obj$cts, na.rm = TRUE)
    if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
    
    return(obj)
    #return(list(cts = include, per = per, clr = clr, bin = bin, nms = colnames(include), meta = obj$meta))
  }
  else {
    if (vrb == TRUE) cat("No OTUs with prevalence <", cutoff, "\b% so no \"Other\" category created.\n")
    return(obj)
  }
}

otu.filter_unc <- function(obj, transform = TRUE, vrb = TRUE) {
  cts = obj$cts
  starting_cts = sum(cts)
  cat("Filter out \"Unclassified\"\n")
  if ("Unclassified" %in% colnames(cts)) {
    cat("Found \"Unclassified\" category in input data\n")
    uncounts = cts[,"Unclassified"]
    if ("Other" %in% colnames(cts)) {
      cat("Found \"Other\" category in input data.\n")
      cts[ ,"Other"] <- cts[ ,"Other"] + uncounts
    } else {
      other = as.matrix(uncounts)
      colnames(other) <- "Other"
      
      obj$orignms <- append(obj$orignms, "Other")
      obj$phynms <- append(obj$phynms, "Other")
      cts <- cbind(cts, other)
      cat("Created new \"Other\" category.\n")
    } 
    
    new_cts = cts[ , colnames(cts) != "Unclassified"]
    obj$orignms = obj$orignms[colnames(cts) != "Unclassified"]
    obj$phynms = obj$phynms[colnames(cts) != "Unclassified"]
    
    per = clr = bin = ""
    if (transform == TRUE) {
      per = per_transform(new_cts)
      clr = clr_transform(new_cts)
      bin = bin_transform(new_cts)
    }
    if (vrb == TRUE) {
      #cat("Collapsed",below_cutoff,"OTUs to \"Other\" otu category.\n")
      cat("Converted",sum(uncounts),"counts to \"Other\" otu category.\n")
      cat("Remaining OTUS:",ncol(new_cts)," (Including \"Other\").\n\n")
    } else {  # No Unclassified
      cat("No \"Unclassified\" category in input data\n")
    }
    
    obj$cts = new_cts
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(new_cts)
    
  } else {
    cat("No \"Unclassified\" category found: data unchanged")
  }
  
  ending_cts = sum(obj$cts)
  if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
  return(obj)
}

otu.move_other <- function(obj, vrb = TRUE) {
  cts = obj$cts
  starting_cts = sum(cts)
  cat("Move \"Other\"\n")
  if ("Other" %in% colnames(cts)) {
    cat("Found \"Other\" category in input data\n")
    othercounts = cts[,"Other"]
    
    cts = cts[ , colnames(cts) != "Other"]
    obj$orignms = obj$orignms[colnames(cts) != "Other"]
    obj$phynms = obj$phynms[colnames(cts) != "Other"]
    obj$orignms <- append(obj$orignms, "Other")
    obj$phynms <- append(obj$phynms, "Other")
    
    other = as.matrix(othercounts)
    colnames(other) <- "Other"
    cts <- cbind(cts, other)
    
    per = clr = bin = ""
    per = per_transform(cts)
    clr = clr_transform(cts)
    bin = bin_transform(cts)
  
    if (vrb == TRUE) {
      #cat("Converted",sum(uncounts),"counts to \"Other\" otu category.\n")
      #cat("Remaining OTUS:",ncol(cts)," (Including \"Other\").\n\n")
    } else {  # No Unclassified
      #cat("No \"Unclassified\" category in input data\n")
    }
    
    obj$cts = cts
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(cts)
  } else {
    cat("No \"Other\" category found: data unchanged")
  }
  
  ending_cts = sum(obj$cts)
  if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
  return(obj)
}

otu.filter_bact <- function(obj, transform = TRUE, vrb = TRUE) {
  cts = obj$cts
  starting_cts = sum(cts)
  cat("Filter out \"Bacteria\"\n")
  if ("Bacteria" %in% colnames(cts)) {
    cat("Found \"Bacteria\" category in input data\n")
    uncounts = cts[,"Bacteria"]
    if ("Other" %in% colnames(cts)) {
      cat("Found \"Other\" category in input data.\n")
      cts[ ,"Other"] <- cts[ ,"Other"] + uncounts
    } else {
      other = as.matrix(uncounts)
      colnames(other) <- "Other"
      
      obj$orignms <- append(obj$orignms, "Other")
      obj$phynms <- append(obj$phynms, "Other")
      cts <- cbind(cts, other)
      cat("Created new \"Other\" category.\n")
    } 
    
    new_cts = cts[ , colnames(cts) != "Bacteria"]
    obj$orignms = obj$orignms[colnames(cts) != "Bacteria"]
    obj$phynms = obj$phynms[colnames(cts) != "Bacteria"]
    
    per = clr = bin = ""
    if (transform == TRUE) {
      per = per_transform(new_cts)
      clr = clr_transform(new_cts)
      bin = bin_transform(new_cts)
    }
    if (vrb == TRUE) {
      #cat("Collapsed",below_cutoff,"OTUs to \"Other\" otu category.\n")
      cat("Converted",sum(uncounts),"counts to \"Other\" otu category.\n")
      cat("Remaining OTUS:",ncol(new_cts)," (Including \"Other\").\n\n")
    } else {  # No Bacteria
      cat("No \"Bacteria\" category in input data\n")
    }
    
    obj$cts = new_cts
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(new_cts)
    
  } else {
    cat("No \"Bacteria\" category found: data unchanged\n")
  }
  
  ending_cts = sum(obj$cts)
  if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
  return(obj)
}

otu.include_otunames <- function(obj, nms, vrb = TRUE, fullnames=TRUE, transform = TRUE) {
  cts = obj$cts

  if (fullnames == TRUE) {
    all_names = obj$orignms
  } else {
    all_names = obj$nms
  }
  
  starting_cts = sum(cts)

  include_v = rep(FALSE, length = ncol(cts))
  for (nm in nms) {
    include = grepl(nm, all_names)
    cat("Found: ", sum(include)," for ",  nm, "\n")
    include_v = include_v | include
  }
  
 # for (i in 1:length(include_v)) {
 #   if (include_v[i] == TRUE) cat(i, ": ", all_names[i], "\n")
 # }
  
  if (sum(include_v) > 0) {
    include = cts[,include_v]
    exclude = cts[,!include_v]
    
    obj$orignms = obj$orignms[include_v]
    obj$phynms = obj$phynms[include_v]
    
    if (class(exclude) == "matrix") other = apply(exclude, MARGIN = 1, sum)
    else if (class(exclude) == "integer") other = exclude
    
    if ("Other" %in% colnames(include)) {
      include[ ,"Other"] <- include[ ,"Other"] + other
      cat("Found \"Other\" category in input data.\n")
    } else {
      other = as.matrix(other)
      colnames(other) <- "Other"
      
      obj$orignms <- append(obj$orignms, "Other")
      obj$phynms <- append(obj$phynms, "Other")
      include <- cbind(include, other)
      cat("Created new \"Other\" category.\n")
    }
    
    if (sum(include, na.rm = TRUE) != sum(obj$cts, na.rm = TRUE)) stop("Error in filtering OTU table!")
    
    per = clr = bin = ""
    if (transform == TRUE) {
      per = per_transform(include)
      clr = clr_transform(include)
      bin = bin_transform(include)
    }
    if (vrb == TRUE) {
#      cat("Retained", sum(include_v),"OTUs.\n")
      cat("Collapsed", sum(!include_v),"OTUs to \"Other\" otu category.\n")
      cat("Converted",sum(exclude),"counts to \"Other\" otu category.\n")
      cat("Remaining OTUS:",ncol(include)," (Including \"Other\").\n\n")
    }
    obj$cts = include
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(include)
    
    ending_cts = sum(obj$cts, na.rm = TRUE)
    if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
    
    return(obj)
  } else {
    if (vrb == TRUE) cat("No OTUs with these names so no \"Other\" category created.\n")
    return(obj)
  }
}

otu.exclude_otunames <- function(obj, nms, vrb = TRUE, fullnames=TRUE, transform = TRUE) {
  cts = obj$cts
  
  if (fullnames == TRUE) {
    all_names = obj$orignms
  } else {
    all_names = obj$nms
  }
  
  starting_cts = sum(cts)
  exclude_v = rep(FALSE, length = ncol(cts))
  for (nm in nms) {
    exclude_nm = grepl(nm, all_names)
    cat("Found: ", sum(exclude_nm)," for ",  nm, "\n")
    exclude_v = exclude_v | exclude_nm
  }
  
  
  if (sum(exclude_v) > 0) {
    include = cts[,!exclude_v]
    exclude = cts[,exclude_v]
    
    obj$orignms = obj$orignms[!exclude_v]
    obj$phynms = obj$phynms[!exclude_v]
    
    if (class(exclude) == "matrix") other = apply(exclude, MARGIN = 1, sum)
    else if (class(exclude) == "integer") other = exclude
    
    if ("Other" %in% colnames(include)) {
      include[ ,"Other"] <- include[ ,"Other"] + other
      cat("Found \"Other\" category in input data.\n")
    } else {
      other = as.matrix(other)
      colnames(other) <- "Other"
      
      obj$orignms <- append(obj$orignms, "Other")
      obj$phynms <- append(obj$phynms, "Other")
      include <- cbind(include, other)
      cat("Created new \"Other\" category.\n")
    }
    
    if (sum(include, na.rm = TRUE) != sum(obj$cts, na.rm = TRUE)) stop("Error in filtering OTU table!")
    
    per = clr = bin = ""
    if (transform == TRUE) {
      per = per_transform(include)
      clr = clr_transform(include)
      bin = bin_transform(include)
    }
    if (vrb == TRUE) {
#      cat("Retained", sum(!exclude_v),"OTUs.\n")
      cat("Collapsed",sum(exclude_v),"OTUs to \"Other\" otu category.\n")
      cat("Converted",sum(exclude),"counts to \"Other\" otu category.\n")
      cat("Remaining OTUS:",ncol(include)," (Including \"Other\").\n\n")
    }
    obj$cts = include
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(include)
    
    ending_cts = sum(obj$cts, na.rm = TRUE)
    if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
    
    return(obj)
  } else {
    if (vrb == TRUE) cat("No OTUs with these names so no \"Other\" category created.\n")
    return(obj)
  }
}

otu.prev <- function(vec) {
  countif = sum(vec > 0)  #obj >0 return vector of TRUE/FALSE (encoded as 1/0) and sum counts the TRUE vals
  count = length(vec)
  prev = (100.0*countif)/count
  return(prev)
}

otu.avg_ra <- function(obj) {
  ras = apply(obj$per, MARGIN = 2, sum) / nrow(obj$per)
  return(ras)
}

otu.filter_ra <- function(obj, cutoff = cutoff, transform = TRUE, vrb = TRUE) {
  if ((cutoff > 100) | (cutoff < 0)) stop("Cutoff must be between 0 and 100 (%)")
  if (cutoff == 0) {
    if (vrb == TRUE) cat("Object unchanged since relative abundance cutoff set to 0%")
    return(obj)
  }
  starting_cts = sum(obj$cts, na.rm = TRUE)
  
  if (vrb == TRUE) cat("Relative abundance cutoff:",cutoff,"\b% (i.e., at least one library must have RA >",cutoff,"\b% to keep OTU).\n")
  cutoff = cutoff / 100.0 #convert from percentage to [0.0, 1.0]
  cts = obj$cts
  per = obj$per
  if (class(per) != "matrix") per = per_transform(cts)
  ras = apply(per, MARGIN = 2, otu.maxabund)
  below_cutoff = sum(ras < cutoff, na.rm = TRUE) #store TRUE/FALSE for whether an OTU has a max abundance < cutoff
  if (below_cutoff == ncol(obj$cts)) stop("All OTUs would be excluded with this cutoff. Halting operation!")
        
  if (sum(below_cutoff, na.rm = TRUE) > 0) {  
    include = cts[,ras >= cutoff]
    exclude = as.matrix(cts[,ras < cutoff])
    above_cutoff = (ras >= cutoff)
    obj$orignms = obj$orignms[ras >= cutoff]
    obj$phynms = obj$phynms[ras >= cutoff]
    if (nrow(exclude) > 1) {
      other = as.matrix(apply(exclude, MARGIN = 1, sum))
    } else if (nrow(exclude) == 1)
      other = exclude
    else {
      cat("Should not arrive here!")
    }
  
    if ("Other" %in% colnames(include)) {
      include[ ,"Other"] <- include[ ,"Other"] + other
      cat("Found \"Other\" category in input data.\n")
    } else {
      other = as.matrix(other)
      colnames(other) <- "Other"
      include <- cbind(include, other)
      obj$orignms <- append(obj$orignms, "Other")
      obj$phynms <- append(obj$phynms, "Other")
      cat("Created new \"Other\" category.\n")
    }

    if (sum(include) != sum(obj$cts, na.rm = TRUE)) stop("Error in filtering OTU table!")
    
    per = clr = bin = ""
    if (transform == TRUE) {
      per = per_transform(include)
      clr = clr_transform(include)
      bin = bin_transform(include)
    }
    if (vrb == TRUE) {
      cat("Collapsed",below_cutoff,"OTUs to \"Other\" otu category.\n")
      cat("Converted",sum(exclude),"counts to \"Other\" otu category.\n")
      cat("Remaining OTUS:",ncol(include)," (Including \"Other\").\n\n")
    }
    obj$cts = include
    obj$per = per
    obj$clr = clr
    obj$bin = bin
    obj$nms = colnames(include)
    ending_cts = sum(obj$cts, na.rm = TRUE)
    if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
    
    return(obj)
    #return(list(cts = include, per = per, clr = clr, bin = bin, nms = colnames(include), meta = obj$meta))
  }
  else {
    if (vrb == TRUE) cat("No OTUs with max. relative abundance <", 100*cutoff, "\b% so no \"Other\" category created.\n")
    return(obj)
  }
}

otu.maxabund <- function(vec) {
  return(max(vec, na.rm = TRUE))
}

otu.plot_counts <- function(obj, mask = "") {
  if (length(mask) <= 1) mask <- !vector(mode = "logical", length = otu.nsubj(obj)) #new vector is filled with "FALSE", so take inverse
  Reads = otu.read_cts(obj)
  Reads = Reads[mask]
  #print(cts)
  par(mfrow=c(1,2))
  
  density = FALSE
  if (density == TRUE) {
    hist(Reads, breaks=20, freq=FALSE, xlim = c(min(Reads),max(Reads)), col = "turquoise", border="blue")
    lines(density(Reads), col="red")
    hist(log10(Reads), breaks=20, freq=FALSE, xlim = c(min(log10(Reads)-0.5),max(log10(Reads)+0.5)), col = "yellow", border="red")
    lines(density(log10(Reads)), col="blue")
  } else {
    hist(Reads, breaks=20, freq=TRUE, xlim = c(min(Reads),max(Reads)), col = "turquoise", border="blue")
    hist(log10(Reads), breaks=20, freq=TRUE, xlim = c(min(log10(Reads)-0.5),max(log10(Reads)+0.5)), col = "yellow", border="red")
  }
  return(summary(Reads))
}

otu.drop_zero_counts <- function(obj, vrb = TRUE) {
  cts = obj$cts
  sums = apply(cts, MARGIN = 2, sum)
  include = cts[,sums > 0]
  #exclude = cts[,sums == 0]
  dropped = ncol(cts) - ncol(include)
  if (dropped > 0) {
    if (dropped >1 & vrb == TRUE) cat("Removed",dropped,"OTUs with zero counts.\n")
    else if (vrb == TRUE) cat("Removed 1 OTU with zero counts.\n")
    if (vrb == TRUE) cat("Remaining OTUS:",ncol(include),"\n")
    
    obj$cts = include
    obj$nms = colnames(include)
    obj$orignms = (obj$orignms)[sums > 0]
    obj$phynms = (obj$phynms)[sums > 0]
    obj = otu.mul_transform(obj)
  } else if (vrb == TRUE) cat("Drop zeros: No OTUs with zero counts found. Return object unchanged.\n")
    
  return(obj)
}

otu.clear_transform <- function(obj = obj) {
  obj$per = ""
  obj$clr = ""
  obj$bin = ""
  return(obj)
}

otu.plot_pvalue_heat <- function(mat, fdr = FALSE, maxp = 0.05, main = "") {
  
  cut = -log10(maxp)
  logmat = mat
  n = ncol(logmat)
  for (i in 1:n) {
    for (j in 1:n) {
      if (mat[i,j] == 0) mat[i,j] = 0.000001
      lp = -log10(mat[i,j])
      if (lp < cut) lp = 0
      logmat[i,j] = lp
    }
  }
  heatmap.2(logmat,col=colorpanel(100,"white","blue","red"),trace="none", dendrogram="none", Rowv=FALSE, Colv = FALSE, main = main)
  
}

otu.mul_transform <- function(obj) {
  obj$per = per_transform(obj$cts)
  obj$clr = clr_transform(obj$cts)
  #print("Here:  otu.mul_transform")
  obj$bin = bin_transform(obj$cts)
  return(obj)
}

otu.shortnames <- function(charvec = charvec) {
  vec = sub("/.../", "/", charvec, ignore.case = FALSE, fixed = TRUE)
  vec = sub("Bacteria/", "", charvec, ignore.case = FALSE, fixed = TRUE)
  
  nms = strsplit(vec, "/", fixed=TRUE)

  phynms = rep("", length(vec))
  otunms = rep("", length(vec))
  for (i in 1:length(nms)) {
    phynm = ""
    otunm = ""
    spl = nms[[i]]
    spllen = length(spl)
    if (spllen > 1) {
      phynm = spl[1]
      phynms[i] = phynm
      phynm = substr(phynm, start = 1, stop = 4)
      otunm = paste(":", spl[spllen], sep = "")
    } else {
      phynm = spl[1]
      phynms[i] = phynm
    }
    otunms[i] = paste(phynm, otunm, sep = "")
  }
  return(list(short = otunms, phynms = phynms))
}


otu.sort_otus <- function(obj, method = "alpha_asc", fullnames = FALSE) {
  
  if (method == "alpha_asc") {
    if (fullnames == TRUE) {
      ord = order(obj$orignms, decreasing = FALSE)
    } else {
      ord = order(colnames(obj$cts), decreasing = FALSE)
    }
  } else if (method == "alpha_desc") {
    if (fullnames == TRUE) {
      ord = order(obj$orignms, decreasing = TRUE)
    } else {
      ord = order(colnames(obj$cts), decreasing = TRUE)
    }
  } else if (method == "ra_asc") {
    abund = otu.avg_ra(obj)
    ord = order(abund, decreasing = FALSE)
  } else if (method == "ra_desc") {
    abund = otu.avg_ra(obj)
    ord = order(abund, decreasing = TRUE)
  }
  
  obj$cts = (obj$cts)[ , ord]
  obj$per = (obj$per)[ , ord]
  obj$clr = (obj$clr)[ , ord]
  obj$bin = (obj$bin)[ , ord]
  obj$nms = (obj$nms)[ord]
  obj$phynms = (obj$phynms)[ord]
  obj$orignms = (obj$orignms)[ord]
  return(obj)
}

otu.topn <- function(obj, n, transform = TRUE, vrb = TRUE) {
  if (n >= otu.notu(obj)) {
    cat("Object unchanged since N >= number of OTUs.\n")
    return(obj)
  }
  if (n <= 0) {
    cat("N must be greater than zero.\n")
    return(obj)
  }
  obj <- otu.sort_otus(obj, method = "ra_desc")
  
  starting_cts = sum(obj$cts, na.rm = TRUE)
  
  cts = obj$cts

  include = as.matrix(cts[,1:n])
  exclude = as.matrix(cts[ , (n+1): ncol(cts)])
  obj$orignms = obj$orignms[1:n]
  obj$phynms = obj$phynms[1:n]
  
  if (nrow(exclude) > 1) {
    other = as.matrix(apply(exclude, MARGIN = 1, sum))
  } else if (nrow(exclude) == 1)
    other = exclude
  else {
    cat("Should not arrive here!")
  }
  
  if ("Other" %in% colnames(include)) {
    include[ ,"Other"] <- include[ ,"Other"] + other
    cat("Found \"Other\" category in input data.\n")
  } else {
    other = as.matrix(other)
    colnames(other) <- "Other"
    include <- cbind(include, other)
    obj$orignms <- append(obj$orignms, "Other")
    obj$phynms <- append(obj$phynms, "Other")
    cat("Created new \"Other\" category.\n")
  }
  
  if (sum(include) != sum(obj$cts, na.rm = TRUE)) stop("Error in filtering OTU table!")
  
  per = clr = bin = ""
  if (transform == TRUE) {
    #cat(ncol(include))
    per = per_transform(include)
    #cat("1")
    clr = clr_transform(include)
    #cat("2")
    bin = bin_transform(include)
    #cat("3")
  }

  if (vrb == TRUE) {
    cat("Collapsed",ncol(exclude),"OTUs to \"Other\" otu category.\n")
    cat("Converted",sum(exclude),"counts to \"Other\" otu category.\n")
    cat("Remaining OTUS:",ncol(include)," (Including \"Other\").\n\n")
  }
  obj$cts = include
  obj$per = per
  obj$clr = clr
  obj$bin = bin
  obj$nms = colnames(include)
  ending_cts = sum(obj$cts, na.rm = TRUE)
  if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
  
  return(obj)
}

otu.plot_manhattan <- function(vals, main = "", sub ="", ylab="", logtrans=TRUE, plot="pval") {
  len = length(vals)
  ind = c(1:len)
  if (plot=="cor") logtrans = FALSE
  if (logtrans == TRUE) y = -log10(vals)
  else y = vals
  if (main == "") main="Manhattan Plot"
  if (ylab == "" & plot=="pval") ylab="-log10(P-value)"
  

  if (plot == "pval") {
    plot(ind, y, pch=4, cex=1, col="blue", main=main, xlab="OTU index", ylab=ylab)
    minor.tick(nx=5, ny=2)
    lines(ind, y, col="blue")
    
    abline(h = -log10(0.01), col="darkgreen", lty=6)
    axis(side=4, at=-log10(0.01), labels="0.01", las=1, cex.axis=0.75, col = "darkgreen")
    abline(h = -log10(0.05), col="darkgreen", lty=6)
    axis(side=4, at=-log10(0.05), labels="0.05", las=1, cex.axis=0.75, col = "darkgreen")
    abline(h = -log10(0.1), col="red", lty=6)
    axis(side=4, at=-log10(0.1), labels="0.1", las=1, cex.axis=0.75, col = "red")
  } else if (plot == "cor") {
    plot(ind, y, pch=4, cex=1, col="blue", main=main, xlab="OTU index", ylab=ylab, ylim=c(-1.0, 1.0))
    minor.tick(nx=5, ny=2)
    lines(ind, y, col="blue")
    abline(h = 0, col="red", lty=6)
    axis(side=4, at=0, labels="0.0", las=1, cex.axis=0.75, col = "red")
    abline(h = 0.5, col="darkgreen", lty=6)
    axis(side=4, at=0.5, labels="0.5", las=1, cex.axis=0.75, col = "darkgreen")
    abline(h = -0.5, col="darkgreen", lty=6)
    axis(side=4, at=-0.5, labels="-0.5", las=1, cex.axis=0.75, col = "darkgreen")
  }
}

otu.plot_pvalue_bar <- function(vals, signs, main = "", sub ="", ylab="", xlab="OTU index", 
                                logtrans=TRUE, plot="pval", lwd = 2, cols = NA) {
  len = length(vals)
  ind = c(1:len)
  if (plot=="cor") logtrans = FALSE
  if (logtrans == TRUE) y = -log10(vals)
  else y = vals
  if (main == "") main="RockyMtn Plot"
  if (ylab == "" & plot=="pval") ylab="-log10(P-value)"
  
  #Convert NAs in -log(p) to 0
  y[is.na(y)] = 0
  signs[is.na(signs)] = TRUE
  
  if (length(signs) > 1) {
    y[signs == "FALSE"] = -1*y[signs == "FALSE"]
  }
  
  #y[is.na(y)] = 0
  if (is.na(cols[1])) {
    cols = c(1:length(y))
    cols[y >= 0] = "blue"
    cols[y < 0] = "red"
  }
  
  if (plot == "pval") {
    
    plot(ind, y, type="h", lwd=lwd, col=cols, main=main, xlab=xlab, ylab=ylab)
    
    minor.tick(nx=5, ny=2)
    
    abline(h = -log10(0.01), col="darkgreen", lty=6)
    axis(side=4, at=-log10(0.01), labels="0.01", las=1, cex.axis=0.75, col = "darkgreen")
    abline(h = -log10(0.05), col="darkgreen", lty=6)
    axis(side=4, at=-log10(0.05), labels="0.05", las=1, cex.axis=0.75, col = "darkgreen")
    abline(h = -log10(0.1), col="red", lty=6)
    axis(side=4, at=-log10(0.1), labels="0.1", las=1, cex.axis=0.75, col = "red")
    
    if (length(signs) > 1) {
      abline(h = -0, col="black", lty=6)
      abline(h = log10(0.01), col="darkgreen", lty=6)
      axis(side=4, at=log10(0.01), labels="0.01", las=1, cex.axis=0.75, col = "darkgreen")
      abline(h = log10(0.05), col="darkgreen", lty=6)
      axis(side=4, at=log10(0.05), labels="0.05", las=1, cex.axis=0.75, col = "darkgreen")
      abline(h = log10(0.1), col="red", lty=6)
      axis(side=4, at=log10(0.1), labels="0.1", las=1, cex.axis=0.75, col = "red")
    }
  } else if (plot == "cor") {
    #plot(ind, y, pch=4, cex=1, col="blue", main=main, xlab="OTU index", ylab=ylab, ylim=c(-1.0, 1.0))
    
    plot(ind, y, type="h", lwd=lwd, col=cols, main=main, xlab=xlab, ylab=ylab, ylim=c(-1.0, 1.0))
    
    minor.tick(nx=5, ny=2)
    #lines(ind, y, col="blue")
    abline(h = 0, col="red", lty=6)
    axis(side=4, at=0, labels="0.0", las=1, cex.axis=0.75, col = "red")
    abline(h = 0.5, col="darkgreen", lty=6)
    axis(side=4, at=0.5, labels="0.5", las=1, cex.axis=0.75, col = "darkgreen")
    abline(h = -0.5, col="darkgreen", lty=6)
    axis(side=4, at=-0.5, labels="-0.5", las=1, cex.axis=0.75, col = "darkgreen")
  }
}

col_pch_assign_factor <- function(in_factor, num.colors = -1) {
  levs = levels(in_factor)
  levs_len = length(levs)
  if (num.colors < levs_len) num.colors = levs_len #num.colors must be >= number of levels
 
  cols = rich.colors(num.colors+1)
  cols <- cols[-1]
  #cols = rainbow(levs_len+1)
  #cols = heat.colors(levs_len+1)
  pchs = rep(c(21, 22, 23, 24, 25), 1 + (levs_len-1)/5)
  szs = rep(c(1.25, 1.25, 1.25, 1, 1), 1 + (levs_len-1)/5)
  len = length(in_factor)
  valvec = in_factor
  colvec = rep("#000000", len)
  pchvec = rep(0, len)
  szvec = rep(0, len)
  for (i in 1:levs_len) {
    colvec[in_factor == levs[i]] <- cols[i]
    pchvec[in_factor == levs[i]] <- pchs[i]
    szvec[in_factor == levs[i]] <- szs[i]
  }
  return(list(valvec = in_factor, colvec = colvec, pchvec = pchvec, cexvec = szvec, lgdcol = cols, lgdpch = pchs, lgdlvl = levs))
}

col_pch_assign_numeric <- function(in_vec, num.colors = -1) {
  cols = colorassign(in_vec, norm = "TRUE")
  colvec = cols$cols
  pchs = pchvec = rep(21, length(in_vec))
  #levs = szvec = rep(1, length(in_vec))
  szvec = 0.2+cols$norm/25
  # cex=0.2+rfcountc$norm[age[maskall] == "02W"]/20
  
  return(list(valvec = in_vec, colvec = colvec, pchvec = pchvec, cexvec = szvec, lgdcol = "", lgdpch = "", lgdlvl = ""))
}


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), lgd.title=lgd.title) {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,5), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=lgd.title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,1,y+1/scale, col=lut[i], border=NA)
  }
}

colorassign <- function(dat, norm = "TRUE") {
  colfunc = colorRampPalette(c("blue", "white", "red"))
  #colfunc = colorRampPalette(c("blue", "lightblue", "red"))
  hc = colfunc(101)

  outlist = rep("#000000", length(dat))
  len = length(dat)
  #print(dat)
  if (norm == "TRUE") {
    normdat = 100*(dat - min(dat, na.rm ="TRUE"))/(max(dat, na.rm ="TRUE") - min(dat, na.rm ="TRUE"))
    dat = normdat
  } 
  #print(dat)
  for (i in 1:len) {
    val = dat[i]
    outlist[i] = hc[val+1]
  }
  return(list(cols = outlist, norm=dat))
}