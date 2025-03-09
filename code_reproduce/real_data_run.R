# This file is for running our method for LncExpDB Brain example
# Run on Zaratan

# Update log:
#   08/24: add HEK and MCF example 
#   08/05: alpha=1e-3 for Brain; also rerun for Smart to finalize. 
#   06/13: for Brain, add prior, replace ncRNA&gene names 
#   06/07: for Brain, set niter to 5,000,000
#   05/11: add Fibro example; set niter to 5,000,000

rm(list=ls())

args <- as.numeric(commandArgs(TRUE))
# boot_id <- args[1]
if (args[1]==1) {eg <- "Brain"}
if (args[1]==2) {eg <- "Fibro"}
alpha <- args[2]

library(MASS) ## ginv() from rpcor()
# library(bnlearn) ## mmpc()
library(pcalg) 
library(graph)
library(BiDAG) ## sampleBN()
library(igraph)
library(CCA)
library(CCP)
library(Rcpp)
library(RcppEigen)

# eg <- "Brain"
# source("~/my_code/causal/helper_fun.R")
source("helper_fun.R") ## local
# sourceCpp("~/my_code/causal/helper_rcpp.cpp")
sourceCpp("helper_rcpp.cpp") ## local 


########## FUNCTIONs ##########

## FUNCTION: Fisher's Z transformation and p-value for correlation
fisher.ztrans.cor <- function(r,n) { ## r is correlation
  z  <-  0.5*(log(1+r)-log(1-r))
  pvalue <- pnorm(-abs(z), mean=0, sd=(1/sqrt(n-3)))*2
  return(pvalue)
}

## FUNCTION: Fisher's Z transformation and p-value for partial correlation
fisher.ztrans <- function(r,n,k) { ## r is partial correlation
  z  <-  0.5*(log(1+r)-log(1-r))
  test.stat <- sqrt(n-k-3)*z
  pvalue <- pnorm(-abs(test.stat))*2
  return(pvalue)
}

## FUNCTION: calculate partial correlation by matrix multiplication
partial.cor <- function(x,y,z) {
  res1 <- cal.res(x,z)
  res2 <- cal.res(y,z)
  pcor <- as.numeric(cor(res1, res2))
  # return(as.numeric(cor(res1,res2)))
  return(list(pcor=pcor, res1=res1, res2=res2))
}


## FUNCTION: calculate residuals in linear regression by matrix multiplication
cal.res <- function(y,X) {
  # A <- solve(t(X) %*% X)
  # A <- svd.inverse(t(X) %*% X)
  A <- ginv(t(X) %*% X)
  beta <- A %*% t(X) %*% y
  res <- (y - X %*% beta)
  if (sd(res)==0) {res <- res + rnorm(length(res), sd=1e-50)}
  return(res)
}

threshold.pvalue <- function(pvalue, threshold=1e-5) {
  pvalue_threshold <- matrix(0,nrow(pvalue),ncol(pvalue))
  pvalue_threshold[(abs(pvalue) < threshold)] <- 1
  return(pvalue_threshold)
}


rpcor <- function(X, Y, alpha=1e-5, alpha0=1e-2, max.reach=5, adj_x, adj_y,
                  nb_mth=1, xcor.thres=0.8, ycor.thres=0.8, d=20) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  sig_rpcor <- NULL
  noi_rpcor <- NULL
  
  cor.X <- cor(X)
  cor.Y <- cor(Y)
  
  nb_lst_x <- vector("list", p)
  nb_lst_y <- vector("list", q)
  for (i in 1:p) { ## get neighbor sets for each X 
    # nb_lst_x[[i]] <- ifelse(nb_mth==1, setdiff(which(adj_x[i,]!=0), i), 
    #                         setdiff(which(abs(cor.X[i,])>xcor.thres), i))
    if (nb_mth==1) {nb_lst_x[[i]] <- setdiff(which(adj_x[i,]!=0), i)}
    if (nb_mth==2) {
      cor_pos_x <- setdiff(which(abs(cor.X[i,])>xcor.thres), i)
      if (length(cor_pos_x)>d) {
        topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
        cor_pos_x <- cor_pos_x[topd_pos]
      }
      nb_lst_x[[i]] <- cor_pos_x
    }
  }
  
  for (j in 1:q) { ## get neighbor sets for each Y
    # nb_lst_y[[j]] <- ifelse(nb_mth==1, setdiff(which(adj_y[j,]!=0), j), 
    #                         setdiff(which(abs(cor.Y[j,])>ycor.thres), j))
    if (nb_mth==1) {nb_lst_y[[j]] <-  setdiff(which(adj_y[j,]!=0), j)}
    if (nb_mth==2) {
      cor_pos_y <- setdiff(which(abs(cor.Y[j,])>ycor.thres), j)
      if (length(cor_pos_y)>d) {
        topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
        cor_pos_y <- cor_pos_y[topd_pos]
      }
      nb_lst_y[[j]] <- cor_pos_y
    }
  }
  
  sis_pvalue <- matrix(NA, p, q)
  for (j in 1:q) {
    cor_y <- as.vector(cor(X, Y[,j]))
    pvalue_y <- fisher.ztrans.cor(cor_y, n)
    sis_pvalue[,j] <- pvalue_y
  }
  output_ind <- threshold.pvalue(sis_pvalue, threshold=alpha0)
  print(paste("marginal screening done, left with ",sum(output_ind)," pairs", sep=""))
  
  
  for (l in 1:max.reach) { ## PC-algorithm from high correlation threshold to low threshold
    # xcor.thres <- ycor.thres <- cor.thres[l]
    # max_ord <- 0
    new_output_ind <- output_ind ## update once every iteration (order)
    
    for (j in 1:q) {
      cor_pos_y <- nb_lst_y[[j]]
      # if (nb_mth==2 & length(cor_pos_y)>d) {
      #   topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
      #   cor_pos_y <- cor_pos_y[topd_pos]
      # }
      
      for (i in 1:p) {
        if (output_ind[i,j] == 1) {
          act_pos_x <- which(output_ind[,j]==1)
          act_pos_y <- which(output_ind[i,]==1)
          # cor_pos_x <- setdiff(which(abs(cor.X[i,])>xcor.thres), i)
          # if (length(cor_pos_x) > d) {
          #   topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
          #   cor_pos_x <- cor_pos_x[topd_pos]
          # }
          cor_pos_x <- nb_lst_x[[i]]
          # if (nb_mth==2 & length(cor_pos_x)>d) {
          #   topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
          #   cor_pos_x <- cor_pos_x[topd_pos]
          # }
          
          nb_pos_x <- intersect(act_pos_x, cor_pos_x)
          nb_pos_y <- intersect(act_pos_y, cor_pos_y)
          
          if ((length(nb_pos_x)+length(nb_pos_y)) > 0) {
            Z <- cbind(X[,nb_pos_x], Y[,nb_pos_y])
            num_nb <- ncol(Z)
            if (num_nb >= l) {
              all_subset_ind <- combn(1:num_nb, l)
              cl <- 1
              one_edge <- 1
              while(cl<=ncol(all_subset_ind) & one_edge==1) {
                sub_Z <- as.matrix(Z[,all_subset_ind[,cl]]) ## condition set
                one_pcor <- partial.cor(X[,i], Y[,j], sub_Z)$pcor
                # one_pcor <- partial_cor(X[,i], Y[,j], sub_Z)
                pvalue <- fisher.ztrans(one_pcor, n, l)
                if (pvalue > alpha) {one_edge <- 0}
                cl <- cl+1
              }
              new_output_ind[i,j] <- one_edge
            } ## end if number of neighbord greater than l
          } ## end if has neighbord
        } ## end if edge active
      } ## end for on 1:p
    } ## end for on 1:q
    
    output_ind <- new_output_ind
    print(paste(l,"-th order screening done, left with ",sum(output_ind)," pairs", sep=""))
  } ## end loop on l
  
  return(output_ind)
  # return(list(ind=output_ind, x=summary(sapply(nb_lst_x, length)), 
  #        y=summary(sapply(nb_lst_y, length)))) ## how many neighbors in average 
}  

rpcor.rcpp <- function(X, Y, alpha=1e-5, alpha0=1e-2, max.reach=5, adj_x, adj_y,
                       nb_mth=1, xcor.thres=0.8, ycor.thres=0.8, d=20) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  sig_rpcor <- NULL
  noi_rpcor <- NULL
  
  cor.X <- cor(X)
  cor.Y <- cor(Y)
  
  nb_lst_x <- vector("list", p)
  nb_lst_y <- vector("list", q)
  for (i in 1:p) { ## get neighbor sets for each X 
    # nb_lst_x[[i]] <- ifelse(nb_mth==1, setdiff(which(adj_x[i,]!=0), i), 
    #                         setdiff(which(abs(cor.X[i,])>xcor.thres), i))
    if (nb_mth==1) {nb_lst_x[[i]] <- setdiff(which(adj_x[i,]!=0), i)}
    if (nb_mth==2) {
      cor_pos_x <- setdiff(which(abs(cor.X[i,])>xcor.thres), i)
      if (length(cor_pos_x)>d) {
        topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
        cor_pos_x <- cor_pos_x[topd_pos]
      }
      nb_lst_x[[i]] <- cor_pos_x
    }
  }
  
  for (j in 1:q) { ## get neighbor sets for each Y
    # nb_lst_y[[j]] <- ifelse(nb_mth==1, setdiff(which(adj_y[j,]!=0), j), 
    #                         setdiff(which(abs(cor.Y[j,])>ycor.thres), j))
    if (nb_mth==1) {nb_lst_y[[j]] <-  setdiff(which(adj_y[j,]!=0), j)}
    if (nb_mth==2) {
      cor_pos_y <- setdiff(which(abs(cor.Y[j,])>ycor.thres), j)
      if (length(cor_pos_y)>d) {
        topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
        cor_pos_y <- cor_pos_y[topd_pos]
      }
      nb_lst_y[[j]] <- cor_pos_y
    }
  }
  
  sis_pvalue <- matrix(NA, p, q)
  for (j in 1:q) {
    cor_y <- as.vector(cor(X, Y[,j]))
    pvalue_y <- fisher.ztrans.cor(cor_y, n)
    sis_pvalue[,j] <- pvalue_y
  }
  output_ind <- threshold.pvalue(sis_pvalue, threshold=alpha0)
  print(paste("marginal screening done, left with ",sum(output_ind)," pairs", sep=""))
  
  
  for (l in 1:max.reach) { ## PC-algorithm from high correlation threshold to low threshold
    # xcor.thres <- ycor.thres <- cor.thres[l]
    # max_ord <- 0
    new_output_ind <- output_ind ## update once every iteration (order)
    
    for (j in 1:q) {
      cor_pos_y <- nb_lst_y[[j]]
      # if (nb_mth==2 & length(cor_pos_y)>d) {
      #   topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
      #   cor_pos_y <- cor_pos_y[topd_pos]
      # }
      
      for (i in 1:p) {
        if (output_ind[i,j] == 1) {
          act_pos_x <- which(output_ind[,j]==1)
          act_pos_y <- which(output_ind[i,]==1)
          # cor_pos_x <- setdiff(which(abs(cor.X[i,])>xcor.thres), i)
          # if (length(cor_pos_x) > d) {
          #   topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
          #   cor_pos_x <- cor_pos_x[topd_pos]
          # }
          cor_pos_x <- nb_lst_x[[i]]
          # if (nb_mth==2 & length(cor_pos_x)>d) {
          #   topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
          #   cor_pos_x <- cor_pos_x[topd_pos]
          # }
          
          nb_pos_x <- intersect(act_pos_x, cor_pos_x)
          nb_pos_y <- intersect(act_pos_y, cor_pos_y)
          
          if ((length(nb_pos_x)+length(nb_pos_y)) > 0) {
            Z <- cbind(X[,nb_pos_x], Y[,nb_pos_y])
            num_nb <- ncol(Z)
            if (num_nb >= l) {
              all_subset_ind <- combn(1:num_nb, l)
              cl <- 1
              one_edge <- 1
              while(cl<=ncol(all_subset_ind) & one_edge==1) {
                sub_Z <- as.matrix(Z[,all_subset_ind[,cl]]) ## condition set
                # one_pcor <- partial.cor(X[,i], Y[,j], sub_Z)$pcor
                one_pcor <- partial_cor(X[,i], Y[,j], sub_Z)
                pvalue <- fisher.ztrans(one_pcor, n, l)
                if (pvalue > alpha) {one_edge <- 0}
                cl <- cl+1
              }
              new_output_ind[i,j] <- one_edge
            } ## end if number of neighbord greater than l
          } ## end if has neighbord
        } ## end if edge active
      } ## end for on 1:p
    } ## end for on 1:q
    
    output_ind <- new_output_ind
    print(paste(l,"-th order screening done, left with ",sum(output_ind)," pairs", sep=""))
  } ## end loop on l
  
  return(output_ind)
  # return(list(ind=output_ind, x=summary(sapply(nb_lst_x, length)), 
  #        y=summary(sapply(nb_lst_y, length)))) ## how many neighbors in average 
}  


rpcor.yy <- function(Y, alpha=1e-5, alpha0=1e-2, max.reach=5,
                     ycor.thres=0.8, d=20) {
  n <- nrow(Y)
  # p <- ncol(X)
  q <- ncol(Y)
  
  
  # cor.X <- cor(X)
  cor.Y <- cor(Y)
  
  # nb_lst_x <- vector("list", p)
  nb_lst_y <- vector("list", q)
  
  for (j in 1:q) { ## get neighbor sets for each Y
    cor_pos_y <- setdiff(which(abs(cor.Y[j,])>ycor.thres), j)
    if (length(cor_pos_y)>d) {
      topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
      cor_pos_y <- cor_pos_y[topd_pos]
    }
    nb_lst_y[[j]] <- cor_pos_y
  }
  
  sis_pvalue <- matrix(NA, q, q)
  for (j in 1:q) {
    # cor_y <- as.vector(cor.Y[,j])
    pvalue_y <- fisher.ztrans.cor(cor.Y[,j], n)
    sis_pvalue[,j] <- pvalue_y
  }
  output_ind <- threshold.pvalue(sis_pvalue, threshold=alpha0)
  diag(output_ind) <- 0
  print(paste("marginal screening done, left with ",sum(output_ind)," pairs", sep=""))
  
  
  for (l in 1:max.reach) { ## PC-algorithm from high correlation threshold to low threshold
    # xcor.thres <- ycor.thres <- cor.thres[l]
    # max_ord <- 0
    new_output_ind <- output_ind ## update once every iteration (order)
    
    for (j in 1:(q-1)) {
      cor_pos_y <- nb_lst_y[[j]]
      
      for (i in (j+1):q) { ## lower triangle 
        
        if (output_ind[i,j] == 1) {
          act_pos_x <- which(output_ind[,j]==1)
          act_pos_y <- which(output_ind[i,]==1)
          
          cor_pos_x <- nb_lst_y[[i]] ## for Y-Y specifically
          
          nb_pos_x <- intersect(act_pos_x, cor_pos_x)
          nb_pos_y <- intersect(act_pos_y, cor_pos_y)
          nb_pos <- setdiff(union(nb_pos_x,nb_pos_y), c(i,j)) ## for Y-Y specifically
          
          if (length(nb_pos) > 0) {
            Z <- as.matrix(Y[,nb_pos]) ## for Y-Y specifically
            num_nb <- ncol(Z)
            if (num_nb >= l) {
              all_subset_ind <- combn(1:num_nb, l)
              cl <- 1
              one_edge <- 1
              while(cl<=ncol(all_subset_ind) & one_edge==1) {
                sub_Z <- as.matrix(Z[,all_subset_ind[,cl]]) ## condition set
                one_pcor <- partial.cor(Y[,i], Y[,j], sub_Z)$pcor
                pvalue <- fisher.ztrans(one_pcor, n, l)
                if (pvalue > alpha) {one_edge <- 0}
                cl <- cl+1
              }
              new_output_ind[i,j] <- new_output_ind[j,i] <- one_edge
            } ## end if number of neighbord greater than l
          } ## end if has neighbord
        } ## end if edge active
      } ## end for on 1:p
    } ## end for on 1:q
    
    output_ind <- new_output_ind
    print(paste(l,"-th order screening done, left with ",sum(output_ind)," pairs", sep=""))
  } ## end loop on l
  
  return(output_ind)
  # return(list(ind=output_ind, x=summary(sapply(nb_lst_x, length)), 
  #        y=summary(sapply(nb_lst_y, length)))) ## how many neighbors in average 
}  

rpcor.yy.rcpp <- function(Y, alpha=1e-5, alpha0=1e-2, max.reach=5,
                          ycor.thres=0.8, d=20) {
  n <- nrow(Y)
  # p <- ncol(X)
  q <- ncol(Y)
  
  
  # cor.X <- cor(X)
  cor.Y <- cor(Y)
  
  # nb_lst_x <- vector("list", p)
  nb_lst_y <- vector("list", q)
  
  for (j in 1:q) { ## get neighbor sets for each Y
    cor_pos_y <- setdiff(which(abs(cor.Y[j,])>ycor.thres), j)
    if (length(cor_pos_y)>d) {
      topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
      cor_pos_y <- cor_pos_y[topd_pos]
    }
    nb_lst_y[[j]] <- cor_pos_y
  }
  
  sis_pvalue <- matrix(NA, q, q)
  for (j in 1:q) {
    # cor_y <- as.vector(cor.Y[,j])
    pvalue_y <- fisher.ztrans.cor(cor.Y[,j], n)
    sis_pvalue[,j] <- pvalue_y
  }
  output_ind <- threshold.pvalue(sis_pvalue, threshold=alpha0)
  diag(output_ind) <- 0
  print(paste("marginal screening done, left with ",sum(output_ind)," pairs", sep=""))
  
  
  for (l in 1:max.reach) { ## PC-algorithm from high correlation threshold to low threshold
    # xcor.thres <- ycor.thres <- cor.thres[l]
    # max_ord <- 0
    new_output_ind <- output_ind ## update once every iteration (order)
    
    for (j in 1:(q-1)) {
      cor_pos_y <- nb_lst_y[[j]]
      
      for (i in (j+1):q) { ## lower triangle 
        
        if (output_ind[i,j] == 1) {
          act_pos_x <- which(output_ind[,j]==1)
          act_pos_y <- which(output_ind[i,]==1)
          
          cor_pos_x <- nb_lst_y[[i]] ## for Y-Y specifically
          
          nb_pos_x <- intersect(act_pos_x, cor_pos_x)
          nb_pos_y <- intersect(act_pos_y, cor_pos_y)
          nb_pos <- setdiff(union(nb_pos_x,nb_pos_y), c(i,j)) ## for Y-Y specifically
          
          if (length(nb_pos) > 0) {
            Z <- as.matrix(Y[,nb_pos]) ## for Y-Y specifically
            num_nb <- ncol(Z)
            if (num_nb >= l) {
              all_subset_ind <- combn(1:num_nb, l)
              cl <- 1
              one_edge <- 1
              while(cl<=ncol(all_subset_ind) & one_edge==1) {
                sub_Z <- as.matrix(Z[,all_subset_ind[,cl]]) ## condition set
                # one_pcor <- partial.cor(Y[,i], Y[,j], sub_Z)$pcor
                one_pcor <- partial_cor(Y[,i], Y[,j], sub_Z)
                pvalue <- fisher.ztrans(one_pcor, n, l)
                if (pvalue > alpha) {one_edge <- 0}
                cl <- cl+1
              }
              new_output_ind[i,j] <- new_output_ind[j,i] <- one_edge
            } ## end if number of neighbord greater than l
          } ## end if has neighbord
        } ## end if edge active
      } ## end for on 1:p
    } ## end for on 1:q
    
    output_ind <- new_output_ind
    print(paste(l,"-th order screening done, left with ",sum(output_ind)," pairs", sep=""))
  } ## end loop on l
  
  return(output_ind)
  # return(list(ind=output_ind, x=summary(sapply(nb_lst_x, length)), 
  #        y=summary(sapply(nb_lst_y, length)))) ## how many neighbors in average 
}  

calc.cc.x <- function(x, y) {
  xx <- as.matrix(x)
  yy <- as.matrix(y)
  cc_val <- cc(xx, yy)$cor
  # pval <- fisher.ztrans.cor(cc_val, n=nrow(y))
  pval <- p.asym(cc_val, nrow(yy), ncol(xx), ncol(yy))$p.value
  return(pval)
}


calc.cc.y <- function(y, x, ym) {
  yy <- as.matrix(y)
  xx <- as.matrix(x)
  yym <- as.matrix(ym)
  
  if (ncol(xx)>0) {
    cc_val1 <- cc(yy, xx)$cor
    pval1 <- p.asym(cc_val1, nrow(yy), ncol(yy), ncol(xx))$p.value
  } else {
    pval1 <- 1
  }
  
  if (ncol(yym)>0) {
    cc_val2 <- cc(yy, yym)$cor
    pval2 <- p.asym(cc_val2, nrow(yy), ncol(yy), ncol(yym))$p.value
  } else {
    pval2 <- 1
  }
  # pval <- fisher.ztrans.cor(cc_val, n=nrow(x))
  return(min(pval1, pval2))
}

cc.scr <- function(X, Y, alpha, ind_xy, ind_yy) {
  p <- ncol(X)
  q <- ncol(Y)
  
  ## screen X node
  node_ind_x <- rep(0, p)
  for (i in 1:p) {
    pos <- which(ind_xy[i,]!=0)
    if (length(pos)>0) {
      one_pval <- calc.cc.x(X[,i], Y[,pos])
      if (one_pval < alpha) {node_ind_x[i] <- 1}
    }
  }
  
  node_ind_y <- rep(0, q)
  for (j in 1:q) {
    # print(j)
    pos1 <- which(ind_xy[,j]!=0)
    pos2 <- which(ind_yy[j,]!=0)
    if ((length(pos1)+length(pos2))>0) {
      one_pval <- calc.cc.y(Y[,j], X[,pos1], Y[,pos2])
      if (one_pval < alpha) {node_ind_y[j] <- 1}
    }
  }
  
  return(list(x=node_ind_x, y=node_ind_y))
}


cc.scr.xy <- function(X, Y, alpha, ind.mat) {
  p <- ncol(X)
  q <- ncol(Y)
  
  ## screen X node
  node_ind_x <- rep(0, p)
  for (i in 1:p) {
    pos <- which(ind.mat[i,]!=0)
    if (length(pos)>0) {
      one_pval <- calc.cc.x(X[,i], Y[,pos])
      if (one_pval < alpha) {node_ind_x[i] <- 1}
    }
  }
  
  ## screen y node
  node_ind_y <- rep(0, q)
  for (j in 1:q) {
    # print(j)
    pos <- which(ind.mat[,j]!=0)
    if (length(pos)>0) {
      one_pval <- calc.cc.x(Y[,j], X[,pos])
      if (one_pval < alpha) {node_ind_y[j] <- 1}
    }
  }
  
  # return(list(x=node_ind_x, y=node_ind_y))
  
  ## remove rows/columns that don't have pairs after node-wise screening
  temp1 <- which(node_ind_x==1)
  temp2 <- which(node_ind_y==1)
  cc_ind_xy <- ind.mat
  cc_ind_xy[which(node_ind_x==0),] <- 0
  cc_ind_xy[,which(node_ind_y==0)] <- 0
  pos_actx <- intersect(temp1, which(rowSums(cc_ind_xy)>0))
  outx <- rep(0, p)
  outx[pos_actx] <- 1
  return(list(x=outx, y=node_ind_y))
}


#################### load data ####################
if (eg=="Brain") {
  # load("~/my_data/LncExpDB_organ_brain_ncrna_gene_fpkm_ppced_0414.RData") ## updated 04/20
  load("LncExpDB_organ_brain_ncrna_gene_fpkm_ppced_0414.RData") ## local
  X <- ncrna_ppced
  Y <- gene_ppced
  
  ## prior info (ind_edb)
  # load("~/my_data/LncExpDB_Brain_lncrna_gene_external_prior_0613.RData")
  load("LncExpDB_Brain_lncrna_gene_external_prior_0613.RData") # local
  # sum(rowSums(ind_edb)>0) # 124
  colnames(X) <- ncname # use LNCI id 
  colnames(Y) <- gname # use gene symbol 
  
  # 06/13: alpha=1e-3
}

if (eg=="Fibro") { ## added 05/04
  # load("~/my_data/Smart_ncrna_gene_count_ppced_0504.RData")
  load("Smart_ncrna_gene_count_ppced_0504.RData") # local 
  X <- ncrna_ppced
  Y <- gene_ppced
  
  X <- X[which(cell_type=="Fibroblasts"),]
  Y <- Y[which(cell_type=="Fibroblasts"),]
  
  # 08/05: alpha=1e-5
}

if (eg=="HEK") { ## added 08/24
  # load("~/my_data/Smart_ncrna_gene_count_ppced_0504.RData")
  load("Smart_ncrna_gene_count_ppced_0504.RData") # local 
  X <- ncrna_ppced
  Y <- gene_ppced
  
  X <- X[which(cell_type=="HEK293T"),]
  Y <- Y[which(cell_type=="HEK293T"),]
  
  # 08/25: alpha=1e-5
}

if (eg=="MCF") { ## added 08/24
  # load("~/my_data/Smart_ncrna_gene_count_ppced_0504.RData")
  load("Smart_ncrna_gene_count_ppced_0504.RData") # local 

  X <- ncrna_ppced[which(cell_type=="MCF7"),]
  Y <- gene_ppced[which(cell_type=="MCF7"),]
  
  ## remove ncRNAs with most zero counts (added 12/06)
  tmp <- apply(X, 2, function(x) mean(x!=0)) > 0.25
  X <- X[, tmp]
  nc_type <- nc_type[tmp]
  
  # 08/25: alpha=1e-3
}

## filter by variance
# top_x <- order(apply(X, 2, var), decreasing=T)[1:floor(0.1*ncol(X))]
# X <- X[,top_x]
# top_y <- order(apply(Y, 2, var), decreasing=T)[1:floor(0.5*ncol(Y))]
# Y <- Y[,top_y]


## center and scale 
X_std <- scale(X)
Y_std <- scale(Y)

## save ncRNA and gene names
# ncname <- colnames(X)
# gname <- colnames(Y)
# save(ncname, gname, 
#      file=paste("~/my_output/causal/", eg, "_early_stage_var_names.RData", sep=""))

#################### hyper parameter ####################
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
d <- round(n/(log(n)))
rho <- 0.3
prior_pen <- 2 ## prior penalty 

# alpha0 <- 1e-2
# alpha <- 1e-6
# alpha2 <- 1e-6
# alpha2 <- 1e-5

print(paste("n =", n))
print(paste("p =", p))
print(paste("q =", q))
print(paste("d =", d))
print(paste("alpha =", alpha))
print(paste("rho =", rho))

set.seed(1234)

#################### Stage 1 ####################

## rPCor for X-Y
t1 <- Sys.time()
## nb_mth=1 is PC, nb_mth=2 is correlation threshold 
# fit_rpcor <- rpcor(X=X_std, Y=Y_std, alpha=alpha, alpha0=alpha, max.reach=5, ## new rPcor
#                      adj_x=adj_xx, adj_y=adj_yy, nb_mth=1, xcor.thres=rho, ycor.thres=rho)
fit_rpcor <- rpcor.rcpp(X=X_std, Y=Y_std, alpha=alpha, alpha0=alpha, max.reach=5,
                   nb_mth=2, xcor.thres=rho, ycor.thres=rho, d=d) ## rPCor
t2 <- Sys.time()
print(t2-t1)

print(sum(rowSums(fit_rpcor)>0))
print(sum(colSums(fit_rpcor)>0))


# save(fit_rpcor,
#      file=paste("~/my_output/causal/", eg, "_early_stage_fit_rpcor_alpha", -log10(alpha), 
#                 "_d", d, "_rho", 10*rho,".RData", sep=""))
# load(paste("~/my_output/causal/", eg, "_early_stage_fit_rpcor_alpha", -log10(alpha), 
#            "_d", d, "_rho", 10*rho,".RData", sep=""))


## node-wise screening: CCA
# capture.output({fit_ccs <- cc.scr(X, Y, alpha2, fit_rpcor, adj_yy)}, file = nullfile())
capture.output({fit_ccs <- cc.scr.xy(X_std, Y_std, alpha, fit_rpcor)}, file = nullfile())
pos_actx <- which(fit_ccs$x==1)
pos_acty <- which(fit_ccs$y==1)
print(paste("left with", sum(fit_ccs$x), "X's"))
print(paste("left with", sum(fit_ccs$y), "Y's"))

## node-wise screening: remove scatter node 
tmp_xy <- fit_rpcor[pos_actx, pos_acty]
p0 <- length(pos_actx)
q0 <- length(pos_acty)
nb_x <- rowSums(tmp_xy)
nb_y <- colSums(tmp_xy)
node_ind_x <- rep(1, p0)
node_ind_y <- rep(1, q0)
for (i in 1:p0) {
  pos <- which(tmp_xy[i,]!=0)
  if((nb_x[i]+max(nb_y[pos]))<=2) {
    node_ind_x[i] <- node_ind_y[pos] <- 0
  }
}
# for (j in 1:q0) {
#   pos <- which(tmp_xy[,j]!=0)
#   if((nb_y[j]+max(nb_x[pos]))<=2) {
#     node_ind_x[pos] <- node_ind_y[j] <- 0
#   }
# }
pos_actx <- pos_actx[node_ind_x==1]
pos_acty <- pos_acty[node_ind_y==1]
print(paste("left with", length(pos_actx), "X's"))
print(paste("left with", length(pos_acty), "Y's"))


X_sub <- X[,pos_actx]
Y_sub <- Y[,pos_acty]
p1 <- length(pos_actx)
q1 <- length(pos_acty)
adj_xy <- fit_rpcor[pos_actx, pos_acty] ## dimensions p1*q1
print(sum(adj_xy))


## PC-skeleton for X-X 
# suff_X <- list(C=cor(X_sub), n=n)
# fit_pc1 <- skeleton(suff_X, indepTest=gaussCItest, alpha=alpha, labels=colnames(X_sub))
# adj_xx <- t(as(fit_pc1, "amat"))
# print(sum(adj_xx))

## PC-skeleton for Y-Y
# t1 <- Sys.time()
# suff_Y <- list(C=cor(Y_sub), n=n)
# fit_pc2 <- skeleton(suff_Y, indepTest=gaussCItest, alpha=alpha, labels=colnames(Y_sub))
# adj_yy2 <- t(as(fit_pc2, "amat"))
# print(sum(adj_yy2))
# t2 <- Sys.time()
# print(t2-t1)

## rPCor for Y-Y
t1 <- Sys.time()
adj_yy <- rpcor.yy.rcpp(Y=Y_sub, alpha=alpha, alpha0=alpha, max.reach=5,
                        ycor.thres=rho, d=d*2) ## rPCor
t2 <- Sys.time()
print(t2-t1)

print(sum(adj_yy)/2)
# print(sum(rowSums(adj_yy)>0))
# print(sum(colSums(adj_yy)>0))

# save(fit_rpcor, adj_xy, adj_yy, pos_actx, pos_acty,
#      file=paste("~/my_output/causal/", eg, "_0613_fit_stage1_alpha", -log10(alpha),
#                 "_d", d, "_rho", 10*rho,".RData", sep=""))
save(fit_rpcor, adj_xy, adj_yy, pos_actx, pos_acty, # local
     file=paste(eg, "_0805_fit_stage1_alpha", -log10(alpha),
                "_d", d, "_rho", 10*rho,".RData", sep=""))

## log: 0614, 0607

#################### Stage 2 ####################
## load results from Stage 1
# load(paste("~/my_output/causal/", eg, "_0421_fit_stage1_alpha", -log10(alpha),
#            "_d", d, "_rho", 10*rho,".RData", sep=""))
X_sub <- X[,pos_actx]
Y_sub <- Y[,pos_acty]
p1 <- length(pos_actx)
q1 <- length(pos_acty)


## build search space for our method (p1*q1)
mat_sp <- matrix(0, p1+q1, p1+q1) 
for (i in 1:p1) {
  ## X-X parents
  # pos_xx <- which(adj_xx[i, ]!=0)
  # if (length(pos_xx)>0) {
  #   mat_sp[pos_xx, i] <- 1
  #   mat_sp[i, pos_xx] <- 1
  # }
  for (j in 1:q1) {
    ## Y-Y parents
    if (i==1) {
      pos_yy <- which(adj_yy[j,]!=0) 
      if (length(pos_yy)>0) {
        mat_sp[p1+pos_yy, p1+j] <- 1
        mat_sp[p1+j, p1+pos_yy] <- 1
      }
    }
    ## X-Y parents
    if (adj_xy[i,j]==1) {mat_sp[i, p1+j] <- 1}  
  }
}
print(sum(mat_sp))

# mat_sp_sub <- mat_sp[c(pos_actx, p+pos_acty), c(pos_actx, p+pos_acty)]
# pos_rmx <- setdiff(1:p, pos_actx)
# pos_rmy <- setdiff(1:q, pos_acty)
# out_s1 <- mat_sp
# out_s1[c(pos_rmx, p+pos_rmy), c(pos_rmx, p+pos_rmy)] <- 0


# # K <- max(colSums(mat_sp))
# # score_data <- scoreparameters("bge", data)
# # mat_bla <- matrix(0,p+q,p+q) ## block list
# # mat_bla[(p+1):(p+q), 1:p] <- 1
# 
K <- max(colSums(mat_sp))
data_sub <- cbind(X_sub, Y_sub)
# score_data_sub <- scoreparameters("bge", data_sub, bgepar=list(edgepf=5))
mat_bla <- matrix(0,p1+q1,p1+q1) ## block list
mat_bla[(p1+1):(p1+q1), 1:p1] <- 1

if (eg=="Brain") {
  ## prepare prior matrix
  tmp <- matrix(prior_pen, p1, q1)
  tmp[ind_edb[pos_actx, pos_acty]==1] <- 1 ## no penalty if in external DB
  mat_prior <- matrix(prior_pen, p1+q1, p1+q1)
  mat_prior[1:p1, (p1+1):(p1+q1)] <- tmp
  
  score_data_sub <- scoreparameters("bge", data_sub, edgepmat=mat_prior)
} else {
  score_data_sub <- scoreparameters("bge", data_sub)
  
}


## Order-MCMC w/ nodewise screening
set.seed(1234)
t1 <- Sys.time()
fit_our <- sampleBN(score_data_sub, algorithm="order", startspace=mat_sp, alpha=alpha,
                    plus1=F, hardlimit=K, iterations=5000000) ## our idea (no plus 1)
t2 <- Sys.time()
print(t2-t1)

dag_our1 <- dag_our <- matrix(0, p+q, p+q)
dag_our[c(pos_actx, p+pos_acty), c(pos_actx, p+pos_acty)] <- as.matrix(getDAG(fit_our))

sum(dag_our)
sum(dag_our[1:p,(p+1):(p+q)]) # level 1
sum(dag_our[(p+1):(p+q), (p+1):(p+q)]) # level 2
sum(rowSums(dag_our[1:p,])>0) # ncRNA node
length(intersect(union(which(rowSums(dag_our)>0), # gene node
                       which(colSums(dag_our)>0)), (p+1):(p+q)))

if (eg %in% c("Fibro","HEK","MCF")) {
  table(nc_type[rowSums(dag_our[1:p,])>0])
  round(table(nc_type[rowSums(dag_our[1:p,])>0])/sum(rowSums(dag_our[1:p,])>0), 3)
}
# sum(rowSums(dag_our[1:p,1:p])>0)
# sum(colSums(dag_our[1:p,(p+1):(p+q)])>0)
# sum(colSums(dag_our[(p+1):(p+q),(p+1):(p+q)])>0)


# save(fit_our, dag_our, fit_rpcor, adj_xy, adj_yy, pos_actx, pos_acty,
#      file=paste("~/my_output/causal/", eg, "_0613_fit_stage2_alpha", -log10(alpha),
#                 "_d", d, "_rho", 10*rho,".RData", sep=""))

# save(fit_our, dag_our, fit_rpcor, adj_xy, adj_yy, pos_actx, pos_acty,
#      file=paste(eg, "_0805_fit_stage2_alpha", -log10(alpha), # local
#                 "_d", d, "_rho", 10*rho,".RData", sep=""))

save(fit_our, dag_our, fit_rpcor, adj_xy, adj_yy, pos_actx, pos_acty,
     file=paste(eg, "_0805_fit_stage2_alpha", -log10(alpha), # local
                "_d", d, "_rho", 10*rho,"_75zero.RData", sep=""))

## log: 0613, 0607

#################### check results ####################
# load(paste("~/my_output/causal/", eg, "_early_stage_fit_stage2_alpha", -log10(alpha),
#            "_d", d, "_rho", 10*rho,".RData", sep=""))
# 
# sum(fit_rpcor)
# sum(rowSums(fit_rpcor)>0)
# sum(colSums(fit_rpcor)>0)
# 
# length(pos_actx)
# length(pos_acty)
# sum(fit_rpcor[pos_actx, pos_acty])
# 
# sum(adj_yy)/2
# length(pos_acty) * (length(pos_acty)-1) / 2
# 
# 
# # sum(dag_our)
# sum(dag_our[1:p,(p+1):(p+q)])
# sum(dag_our[(p+1):(p+q), (p+1):(p+q)])
# 
# # sum(rowSums(dag_our[(p+1):(p+q),])>0 | colSums(dag_our[(p+1):(p+q),(p+1):(p+q)])>0)
# 
# sum(rowSums(dag_our[1:p,])>0)
# length(intersect(union(which(rowSums(dag_our)>0), which(colSums(dag_our)>0)), (p+1):(p+q)))

