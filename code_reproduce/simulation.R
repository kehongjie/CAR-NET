## Update log:
#### 05/21: adjust number of "negative"; fix AGES adj output
#### 06/01: increase number of edges (20*5+50); modify k_pen to (20,15,10,5,1); 
####        modify d=n/log(n)/2; modify rho=0.7
#### 06/02: modify k_pen to (10,8,6,4,2); 
####        modify alpha to (0.0001,0.0005,0.001,0.005,0.01);
####        modify edges (10*10+50)
#### 06/02_v2: switch back to edges (10*5+50);  modify k_pen to (7,6,5,4,3)


bash_arg <- commandArgs(TRUE) #(simu_id, sce, n, ind_tune)
# bash_arg <- c(1, 1, 200, 1) 
# bash_arg <- c(1, 2, 300, 1)

library(bnlearn) ## mmpc()
library(pcalg) 
library(graph)
library(BiDAG) ## sampleBN()
library(igraph)
library(CCA)
library(CCP)

# rm(list=ls())
# set.seed(12345)
# load("~/my_output/causal/all_seed.RData") ## all_seed
load("all_seed.RData") ## all_seed


########## FUNCTIONs ##########
tp.in.rep <- function(x, level) {
  if (level==1) {
    return(sum(x[pos_l1]))
  }
  if (level==2) {
    return(sum(x[pos_l2]))
  }
}

tpr.in.rep <- function(x, level) {
  if (level==1) {
    return(sum(x[pos_l1]) / sum(A[pos_l1]))
  }
  if (level==2) {
    return(sum(x[pos_l2]) / sum(A[pos_l2]))
  }
  if (level==3) {
    return( (sum(x[pos_l1])+sum(x[pos_l2])) / (sum(A[pos_l1])+sum(A[pos_l2])) )
  }
}

fp.in.rep <- function(x, level) {
  if (level==1) {
    return(sum(x[pos_neg_l1]))
  }
  if (level==2) {
    return(sum(x[pos_neg_l2]))
  }
}

ppv.in.rep <- function(x, level) {
  if (level==1) {
    tp <- sum(x[pos_l1])
    ms <- sum(x[union(pos_l1, pos_neg_l1)])
  }
  if (level==2) {
    tp <- sum(x[pos_l2])
    ms <- sum(x[union(pos_l2, pos_neg_l2)])
  }
  if (level==3) {
    tp <- sum(x[pos_l1])+sum(x[pos_l2])
    ms <- sum(x[union(pos_l1, pos_neg_l1)]) + sum(x[union(pos_l2, pos_neg_l2)])
  }
  return(tp/ms)
}


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
  A <- solve(t(X) %*% X)
  beta <- A %*% t(X) %*% y
  res <- (y - X %*% beta)
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


########## hyper-parameter ##########

beta <- 0.5
s2 <- 0.2 ## sigma square
B <- 20
model <- 3.3
all_adj <- vector("list", B)

simu_id <- as.numeric(bash_arg[1])
sce <- ifelse(as.numeric(bash_arg[2])==1, "lowd", "highd")
# tune_all <- matrix(c(0.001,0.005,0.01,0.05,0.1, ## all tunning parameter choice
#                      5,4,3,2,1), ncol=2)
if (sce=="lowd") {
  tune_all <- matrix(c(0.0001,0.0005,0.001,0.005,0.01, ## all tunning parameter choice
                       7,6,5,4,3), ncol=2) ## updated 06/02
} else {
  # tune_all <- matrix(c(1e-5,5e-4,1e-4,5e-3,1e-3), ncol=1) ## all tunning parameter choice
  tune_all <- matrix(c(1e-5,5e-5,1e-4,5e-4,1e-3), ncol=1) ## all tunning parameter choice
}
n <- as.numeric(bash_arg[3])
ind_tune <- as.numeric(bash_arg[4])

alpha  <- tune_all[ind_tune,1]
if (sce=="lowd") {k_pen  <- tune_all[ind_tune,2]}
num_mth <- ifelse(sce=="lowd", 6, 3)
# alpha2 <- as.numeric(bash_arg[5])

d <- round(0.5*n/(log(n))) ## updated 06/01
rho <- 0.7 ## updated 06/01

print(paste("Simulation No.", simu_id))
print(paste("Scenario is ", sce, sep=""))
print(paste("n=", n, sep=""))
print(paste("ind_tune=", ind_tune, sep=""))
print(paste("alpha=", alpha, sep=""))
if (sce=="lowd") {print(paste("k_pen=", k_pen, sep=""))}
print(paste("beta=", beta, sep=""))


if (sce=="lowd") {
  p <- 100
  q <- 100
  # ss <- 0.2
  # ss_xx <- 0.01 ## X-X sparsity
  # ss_yy <- 0.01 ## Y-Y sparsity (updated 06/01)
  # ss_xy <- 0.5 ## X-Y sparsity
  s_rep <- 10 ## number of repetitive structure (updated 06/01)
  s_hub <- 5 ## number of genes that one hub ncRNA regulates, 5*10 or 10*10?
  s_yy <- 50 ## number of gene-gene edges (updated 06/01)
} else {
  p <- 1000
  q <- 1000
  # ss <- 0.2
  # ss_xx <- 0.0001 ## X-X sparsity
  # ss_yy <- 0.02/(10^2) ## Y-Y sparsity
  # ss_xy <- 0.005 ## X-Y sparsity
  s_rep <- 10 ## number of repetitive structure 
  s_hub <- 10 ## number of genes that a hub ncRNA regulates
  s_yy <- 50 ## number of gene-gene edges (updated 06/01)
}


# all_adj <- lapply(vector("list",num_mth), function(x) matrix(0,p+q,p+q))

mat_shd <- matrix(NA, B, num_mth)
tp_l1 <- fp_l1 <- tp_l2 <- fp_l2 <- matrix(NA, B, num_mth)

##############################
## !!!!! SIMULATION starts here
##############################
# for (b in 1:B) {
b <- simu_id
set.seed(all_seed[simu_id])

rep_adj <- vector("list", num_mth)
A <- matrix(0, p+q, p+q) ## adj matrix, row to column
colnames(A) <- rownames(A) <- c(paste("X",1:p,sep=""), paste("Y",1:q,sep=""))


########## generate A matrix  ##########
## III: repetitive hub ncRNAs
if (model %in% c(3.1,3.2,3.3)) {
  ## hub ncRNA structure 
  pos1 <- sort(sample(1:p, s_rep))
  pos2 <- sort(sample(1:q, s_hub*s_rep))
  for (i in 1:s_rep) {
    # A[i, ((i-1)*s_hub+1):((i-1)*s_hub+s_hub)]
    for (j in 1:s_hub) {
      A[pos1[i], p+pos2[(i-1)*s_hub+j]] <- 1
    }
  }
  
  ## Y-Y edges
  if (model == (3.1|3.2)) {
    idx_y <- sample(1:choose(q,2), size=s_yy)
    ct <- 0
    for (i in 1:(q-1)) {
      for (j in (i+1):q) {
        ct <- ct+1
        if (ct %in% idx_y) {A[p+i,p+j] <- 1}
      }
    }
  }
  
  ## Y-Y edges (updated 11/12)
  if (model == 3.3) {
    tmpy <- which(colSums(A[,(p+1):(p+q)])>0)
    idx_y <- sample(1:choose(length(tmpy),2), size=s_yy)
    subA <- matrix(0, length(tmpy), length(tmpy))
    ct <- 0
    for (i in 1:(length(tmpy)-1)) {
      for (j in (i+1):length(tmpy)) {
        ct <- ct+1
        if (ct %in% idx_y) {subA[i,j] <- 1}
      }
    }
    A[p+tmpy, p+tmpy] <- subA
  }
  
  
  ## X-X edges
  if (model == 3.1) {
    idx_x <- sample(1:choose(p,2), size=floor(choose(p,2)*ss_xx))
    ct <- 0
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        ct <- ct+1
        if (ct %in% idx_x) {A[i,j] <- 1}
      }
    }
  }
  
  
} ## end model 3.1/3.2/3.3


########## generate data based on A matrix   ##########
if (model != 2.1) {
  data <- matrix(0, n, p+q)
  for (i in 1:(p+q)) {
    data[,i] <- rnorm(n, sd=sqrt(s2)) + data %*% matrix(beta*A[,i], ncol=1)
  }
}

colnames(data) <- c(paste("X",1:p,sep=""), paste("Y",1:q,sep=""))
data <- as.data.frame(data)
X <- data[,1:p]
Y <- data[,(p+1):(p+q)]

# plot(as(A, "graphNEL"))

## Positive & Negative position
mat_xx <- mat_yy <- mat_xy <- mat_yx <- matrix(0, p+q, p+q)
mat_xx[1:p, 1:p] <- 1
for (kk in 1:p) {mat_xx[kk,kk] <- 0} ## updated 05/21
mat_yy[(p+1):(p+q), (p+1):(p+q)] <- 1
for (kk in (p+1):(p+q)) {mat_yy[kk,kk] <- 0} ## updated 05/21
mat_xy[1:p, (p+1):(p+q)] <- 1
mat_yx[(p+1):(p+q), 1:p] <- 1
pos_l1 <- which(A*mat_xy==1)
pos_l2 <- which(A*(mat_xx+mat_yy)==1)
# pos_neg_l1 <- which(A==0 & mat_xy==1)
pos_neg_l1 <- which(A==0 & (mat_xy+mat_yx)>0) # or mat_xy==1? updated 05/21
# pos_neg_l2 <- which(A==0 & mat_xy!=1)
pos_neg_l2 <- which(A==0 & mat_yy==1) ## updated 05/21
# sum(c(pos_l1, pos_l2, union(pos_neg_l2, pos_neg_l1))) ## should be sum(1:length(A))
pos_tx <- which(rowSums(A[1:p,])>0)
pos_ty <- union(which(colSums(A[,(p+1):(p+q)])>0), which(rowSums(A[(p+1):(p+q),])>0))

print(paste("# pos_l1 =", length(pos_l1)))
print(paste("# pos_l2 =", length(pos_l2)))
print(paste("# pos_tx =", length(pos_tx)))
print(paste("# pos_ty =", length(pos_ty)))

###################################
##### below are running algorithms
###################################


#################### Step 1 #################### 
## rPCor for X-Y
t1 <- Sys.time()
## nb_mth=1 is PC, nb_mth=2 is correlation threshold 
# fit_rpcor <- rpcor(X=X, Y=Y, alpha=alpha, alpha0=alpha, max.reach=5, ## new rPcor
#                      adj_x=adj_xx, adj_y=adj_yy, nb_mth=1, xcor.thres=rho, ycor.thres=rho)
fit_rpcor <- rpcor(X=X, Y=Y, alpha=alpha, alpha0=alpha, max.reach=5,
                   nb_mth=2, xcor.thres=rho, ycor.thres=rho, d=d) ## rPCor
t2 <- Sys.time()
print(t2-t1)
print(sum(rowSums(fit_rpcor)>0))
print(sum(colSums(fit_rpcor)>0))


## node-wise screening: CCA
# capture.output({fit_ccs <- cc.scr(X, Y, alpha2, fit_rpcor, adj_yy)}, file = nullfile())
capture.output({fit_ccs <- cc.scr.xy(X, Y, alpha, fit_rpcor)}, file = nullfile())
pos_actx <- which(fit_ccs$x==1)
pos_acty <- which(fit_ccs$y==1)
print(paste("left with", sum(fit_ccs$x), "X aftre CCA"))
print(paste("left with", sum(fit_ccs$y), "Y after CCA"))

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

pos_actx <- pos_actx[node_ind_x==1]
pos_acty <- pos_acty[node_ind_y==1]
print(paste("left with", length(pos_actx), "X after removing scatters"))
print(paste("left with", length(pos_acty), "Y after removing scatters"))


X_sub <- X[,pos_actx]
Y_sub <- Y[,pos_acty]
p1 <- length(pos_actx)
q1 <- length(pos_acty)
adj_xy <- fit_rpcor[pos_actx, pos_acty] ## dimensions p1*q1
print(paste("Number of X-Y edges is", sum(adj_xy)))


## rPCor for Y-Y
t1 <- Sys.time()
adj_yy <- rpcor.yy(Y=Y_sub, alpha=alpha, alpha0=alpha, max.reach=5,
                   ycor.thres=rho, d=d*2) ## rPCor
t2 <- Sys.time()
print(t2-t1)
print(paste("Number of Y-Y edges is", sum(adj_yy)))
# sum(rowSums(adj_yy)>0)

X_sub <- X[,pos_actx]
Y_sub <- Y[,pos_acty]
p1 <- length(pos_actx)
q1 <- length(pos_acty)

## build search space for our method
mat_sp <- matrix(0, p1+q1, p1+q1) 
for (i in 1:p1) {
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
# sum((mat_sp_sub + A[c(pos_actx, p+pos_acty), c(pos_actx, p+pos_acty)])==2) ## quick check

#################### Step 2 #################### 
mat_bla <- matrix(0,p1+q1,p1+q1) ## block list
mat_bla[(p1+1):(p1+q1), 1:p1] <- 1

K <- max(colSums(mat_sp))
data_sub <- cbind(X_sub, Y_sub)
score_data_sub <- scoreparameters("bge", data_sub)
# score_data_sub <- scoreparameters("bge", data_sub, bgepar=list(edgepf=5))

# ## w/o nodewise screening
# fit_our1 <- sampleBN(score_data, algorithm="order", startspace=mat_sp, plus1=T,
#                     blacklist=mat_bla, hardlimit=K) ## our idea (plus 1)
# fit_our <- sampleBN(score_data, algorithm="order", startspace=mat_sp, alpha=alpha, 
#                      plus1=F, hardlimit=K) ## our idea (no plus 1)

## w/ nodewise screening
t1 <- Sys.time()
fit_our <- sampleBN(score_data_sub, algorithm="order", startspace=mat_sp, alpha=alpha, 
                    plus1=F, hardlimit=K) ## our idea (no plus 1)
t2 <- Sys.time()
print(t2-t1)

# fit_our1 <- sampleBN(score_data_sub, algorithm="order", startspace=mat_sp_sub, alpha=alph, 
#                      plus1=T, blacklist=mat_bla_sub, hardlimit=K1) ## our idea (plus 1)
# t3 <- Sys.time()
# print(t3-t2)

dag_our1 <- dag_our <- matrix(0, p+q, p+q)
dag_our[c(pos_actx, p+pos_acty), c(pos_actx, p+pos_acty)] <- as.matrix(getDAG(fit_our))
# dag_our1[c(pos_actx, p+pos_acty), c(pos_actx, p+pos_acty)] <- as.matrix(getDAG(fit_our1))

# print(compareDAGs(dag_our, A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
print("Results for our method:")
print(c(tpr.in.rep(dag_our, 1), tpr.in.rep(dag_our, 2), 
      ppv.in.rep(dag_our, 1), ppv.in.rep(dag_our, 2)))


#################### hybrid MCMC #################### 
if (sce=="lowd") {
  # mat_om <- matrix(1, p+q, p+q)  ## original search space
  # fit_om <- sampleBN(score_data, algorithm="order",  startspace=mat_om, hardlimit=(p+q)) ## regular orderMCMC
  # print(compareDAGs(getDAG(fit_om), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
  # 
  
  score_data <- scoreparameters("bge", data)
  t1 <- Sys.time()
  fit_im <- sampleBN(score_data, algorithm="order", alpha=alpha, 
                     plus1=F) ## iterative orderMCMC (no plus 1)
  t2 <- Sys.time()
  print(t2-t1)
  
  # fit_im1 <- sampleBN(score_data, algorithm="order", alpha=alpha, 
  #                     plus1=T) ## iterative orderMCMC (plus 1)
  
  # print(compareDAGs(getDAG(fit_im), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
  print("Results for hybrid MCMC:")
  print(c(tpr.in.rep(getDAG(fit_im), 1), tpr.in.rep(getDAG(fit_im), 2), 
          ppv.in.rep(getDAG(fit_im), 1), ppv.in.rep(getDAG(fit_im), 2)))

}


########## PC ##########
suffstat <- list(C=cor(data), n=n)
# fit_pcs <- skeleton(suffstat, indepTest=gaussCItest, alpha=alpha, labels=colnames(data))
# tmp <- t(as(fit_pcs, "amat"))

t1 <- Sys.time()
fit_pc <- pc(suffstat, indepTest=gaussCItest, alpha=alpha, labels=colnames(data))
t2 <- Sys.time()
print(t2-t1)

# print(compareDAGs(t(as(fit_pc, "amat")), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
print("Results for PC:")
print(c(tpr.in.rep(t(as(fit_pc, "amat")), 1), tpr.in.rep(t(as(fit_pc, "amat")), 2), 
        ppv.in.rep(t(as(fit_pc, "amat")), 1), ppv.in.rep(t(as(fit_pc, "amat")), 2)))


########## MMPC ##########
## constraint based: mmpc
t1 <- Sys.time()
fit_mmpc <- mmpc(data, alpha=alpha, undirected=F)
t2 <- Sys.time()
print(t2-t1)

# print(compareDAGs(amat(fit_mmpc), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
print("Results for MMPC:")
print(c(tpr.in.rep(amat(fit_mmpc), 1), tpr.in.rep(amat(fit_mmpc), 2), 
        ppv.in.rep(amat(fit_mmpc), 1), ppv.in.rep(amat(fit_mmpc), 2)))


########## AGES ##########
# fit_ages0 <- ages(data, lambda_min=k_pen)
# fit_ages <- as(as(fit_ages0$essgraph, "graphNEL"), "Matrix")

# library(huge)
# huge.path <- huge(gmG8$x, verbose = FALSE)
# huge.fit <- huge.select(huge.path, verbose = FALSE)
# cig <- as.matrix(huge.fit$refit)
# ges.fit <- ges(score, fixedGaps = !cig, adaptive = "vstructures")

score <- new("GaussL0penObsScore", data) 
skel.fit <- skeleton(suffStat = suffstat, indepTest = gaussCItest, 
                     alpha = alpha, p = ncol(data)) 
skel <- as(skel.fit@graph, "matrix")
fit_ages0 <- ges(score, fixedGaps=!skel, adaptive="triples")
fit_ages <- as(as(fit_ages0$essgraph, "graphNEL"), "Matrix")


# print(compareDAGs(as.matrix(fit_ages), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
print("Results AGES:")
print(c(tpr.in.rep(as.matrix(fit_ages), 1), tpr.in.rep(as.matrix(fit_ages), 2), 
        ppv.in.rep(as.matrix(fit_ages), 1), ppv.in.rep(as.matrix(fit_ages), 2)))


########## HC ##########
if (sce=="lowd") {
  ## HC, score based
  t1 <- Sys.time()
  fit_hc <- hc(data, k=k_pen)
  t2 <- Sys.time()
  print(t2-t1)
  # graphviz.plot(fit_hc, highlight=hlnode)
  # rep_adj[[3]] <- amat(fit_hc)
  # print(compareDAGs(amat(fit_hc), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
  print("Results for HC:")
  print(c(tpr.in.rep(amat(fit_hc), 1), tpr.in.rep(amat(fit_hc), 2), 
          ppv.in.rep(amat(fit_hc), 1), ppv.in.rep(amat(fit_hc), 2)))
  # tpr.in.rep(amat(fit_hc), 1)
  # ppv.in.rep(amat(fit_hc), 1)
  
  # fit_tabu <- tabu(data)
  # # graphviz.plot(fit_tabu, highlight=hlnode)
  # # rep_adj[[4]] <- amat(fit_tabu)
  # print(compareDAGs(amat(fit_tabu), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
  # 
  # # 
  ## hybrid methods
  # fit_mmhc <- mmhc(data)
  # # graphviz.plot(fit_mmhc, highlight=hlnode)
  # # rep_adj[[5]] <- amat(fit_mmhc)
  # print(compareDAGs(amat(fit_mmhc), A, cpdag = TRUE)[c("TPR", "FDR", "SHD")])
}

########## Genie 3 ##########



########## calculate TP&FP ##########
if (sce=="lowd") {
  rep_adj[[1]] <- dag_our
  # rep_adj[[2]] <- dag_our1
  rep_adj[[2]] <- getDAG(fit_im)
  # rep_adj[[4]] <- getDAG(fit_im1)
  rep_adj[[3]] <- t(as(fit_pc, "amat"))
  rep_adj[[4]] <- amat(fit_mmpc)
  rep_adj[[5]] <- amat(fit_hc)
  # rep_adj[[6]] <- t(as.matrix(fit_ages))
  rep_adj[[6]] <- as.matrix(fit_ages) ## updated 05/21
} else {
  rep_adj[[1]] <- dag_our
  rep_adj[[2]] <- t(as(fit_pc, "amat"))
  rep_adj[[3]] <- amat(fit_mmpc)
}



# tp_l1[b,] <- sapply(rep_adj, tp.in.rep, level=1)
# tp_l2[b,] <- sapply(rep_adj, tp.in.rep, level=2)
# fp_l1[b,] <- sapply(rep_adj, fp.in.rep, level=1)
# fp_l2[b,] <- sapply(rep_adj, fp.in.rep, level=2)
# 
# print(tp_l1[b,])
# print(tp_l2[b,])
# print(fp_l1[b,])
# print(fp_l2[b,])

# mat_shd[b,] <- c(compareDAGs(dag_our, A, cpdag = TRUE)["SHD"], 
#                  compareDAGs(dag_our1, A, cpdag = TRUE)["SHD"],
#                  compareDAGs(getDAG(fit_im), A, cpdag = TRUE)["SHD"], 
#                  compareDAGs(getDAG(fit_im1), A, cpdag = TRUE)["SHD"],
#                  compareDAGs(t(as(fit_pc, "amat")), A, cpdag = TRUE)["SHD"])

# all_adj[[b]] <- rep_adj
save(rep_adj, A, pos_l1, pos_l2, pos_neg_l1, pos_neg_l2, pos_tx, pos_ty,
     file=paste("~/my_output/causal/simulation/simu0602v2_b", simu_id, "_", sce, "_n", n, "_pq", p, 
                "_tune", ind_tune,".RData", sep="")) ## on Zaratan


# } ## end of replication


########## visualize and summarize ##########
# apply(tp_l1, 2, mean)
# apply(tp_l2, 2, mean)
# apply(fp_l1, 2, mean)
# apply(fp_l2, 2, mean)
# apply(mat_shd, 2, mean)


# save(tp_l1, tp_l2, fp_l1, fp_l2, 
#      file=paste("~/my_output/causal/simu", simu_id, "_", sce, "_n", n, "_pq", p, 
#                 "_tune", ind_tune,".RData", sep=""))


# ## graph comparison 
# plotdiffs(matrix(as.vector(t(as(fit_pc, "amat"))), nrow=(p+q)), as(A, "graphNEL")) ## PC
# plotdiffs(graphviz.plot(fit_pc2), as(A, "graphNEL"))
# plotdiffs(graphviz.plot(fit_cpc), as(A, "graphNEL"))
# plotdiffs(graphviz.plot(fit_tpc), as(A, "graphNEL"))
# 
# plotdiffs(as(getDAG(fit_omc), "graphNEL"), as(A, "graphNEL")) ## order-MCMC
# plotdiffs(as(getDAG(fit_omc2), "graphNEL"), as(A, "graphNEL"))
# plotdiffs(as(getDAG(fit_omc3), "graphNEL"), as(A, "graphNEL")) ## constraint + order-MCMC
# plotdiffs(as(getDAG(fit_omc4), "graphNEL"), as(A, "graphNEL"))
# 
# plotdiffs(as(getDAG(fit_comc), "graphNEL"), as(A, "graphNEL"))
# 
# # 
# # plotdiffs(as(getDAG(fit_pmc), "graphNEL"), as(A, "graphNEL"))
# 
# plot(fit_pc)
# iplotPC(fit_pc)
# graphviz.plot(fit_mmpc, highlight=hlnode)
# graphviz.plot(fit_hc, highlight=hlnode)
# graphviz.plot(fit_tabu, highlight=hlnode)
# graphviz.plot(fit_mmhc, highlight=hlnode)
# plot(as(getDAG(fit_homc), "graphNEL"))
# plot(as(getDAG(fit_hpar), "graphNEL"))
# 


