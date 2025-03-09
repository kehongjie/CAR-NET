
## Functions needed to run screening + order-MCMC
## Last updated: 06/08/2024

## Update log
####  06/08: add Rcpp version of partial.cor()


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

