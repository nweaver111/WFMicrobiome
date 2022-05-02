#************************
#
#
# v.0.1.0  
# daniel.frank@ucdenver.edu
# 24 April 2017
#
#
#************************


mul_transform <- function(cts) {
  per = per_transform(cts)
  clr = clr_transform(cts)
  bin = bin_transform(cts)
  nms = colnames(cts)
  return(list(cts=cts, per=per, clr=clr, bin=bin, nms=nms))
}

#perform CLR transformation for microbiota data
# adapted from code written by Brandie Wagner 04-03-2013
 
clr_transform =  function(in_mat)
{
  #calculate number of samples
  n_samp <- nrow(in_mat)

  #calculate number of taxa 
  n_org <- ncol(in_mat)

  #define matrices
  ra <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  l_ra <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  clr <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  gm <- rep(0,n_samp)
  c <- rep(0,n_samp)
  total <- rep(0,n_samp)
  
  for (i in 1:n_samp) {
	  #calc weight
    total[i] = sum(in_mat[i,])
	  c[i] <- 1/total[i]

  	for (j in 1:n_org) {
	    #transform counts to relative abundance and add small amount and log transform
	    ra[i,j] <- (in_mat[i,j]+c[i])/(total[i]+c[i]*n_org)
      #ra[i,j] <- (in_mat[i,j]+c[i])/(total[i]+c[i]*n_org)
	    l_ra[i,j] <- log(ra[i,j])
	  }
	  #calculate gm
	  gm[i]=mean(l_ra[ i,])
    clr[ i,] <- l_ra[i ,] - gm[i]
  }
  
  colnames(clr) = colnames(in_mat)
  rownames(clr) = rownames(in_mat)
  
  return(clr)
}

#perform relative abundance transformation for microbiota data
# code written by Dan Frank

per_transform =  function(in_mat)
{
  #calculate number of samples
  n_samp <- nrow(in_mat)

  #calculate number of taxa 
  n_org <- ncol(in_mat)

  #define matrices
  ra <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  l_ra <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  clr <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  gm <- rep(0,n_samp)
  c <- rep(0,n_samp)
  total <- rep(0,n_samp)
  
  for (i in 1:n_samp) {
	  #calc weight
    total[i] = sum(in_mat[i,], na.rm = TRUE)

  	for (j in 1:n_org) {
	    ra[i,j] <- in_mat[i,j]/total[i]
	  }
  }
  
  colnames(ra) = colnames(in_mat)
  rownames(ra) = rownames(in_mat)
  return(ra)
}

#perform binary transformation for microbiota data
# code written by Dan Frank

bin_transform =  function(in_mat)
{
  #calculate number of samples
  n_samp <- nrow(in_mat)
  #calculate number of taxa 
  n_org <- ncol(in_mat)

  #define matrices
  bin <- matrix (rep(0,n_org*n_samp), n_samp, n_org)
  
  for (i in 1:nrow(in_mat)) {
  	for (j in 1:ncol(in_mat)) {
      if (in_mat[i,j] > 0) bin[i,j] = 1
	    else bin[i,j] = 0
	  }
  }
  
  colnames(bin) = colnames(in_mat)
  rownames(bin) = rownames(in_mat)
  return(bin)
}
