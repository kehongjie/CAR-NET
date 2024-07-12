library(EnvStats)
library(stringr)

#####reference databases
gene_transcript_0 = read.csv("./data/reference/ENSG_lnci.csv")[,-1]
HSAL_ENSG_lnci = read.csv("./data/reference/HSAL_lnci.csv")[,-1]
HGNC_ENSG_lnci = read.csv("./data/reference/HGNC_lnci.csv")[,-1]

##########internal map functions:
#function: ENSG/ENST(without extension) to lncipedia ID
map_ID_1_0 = function(x){
  if(!is.character(x)){
    return(print("Error: should enter a string"))
  }
  x_0 = word(x, sep = "\\.")
  if(x_0 %in% gene_transcript_0[,2]){
    ind_x_0 = which(gene_transcript_0[,2] == x_0)
    return(gene_transcript_0[,1][ind_x_0])
  }else{
    return(NA)
  }
}
#function: HSALN lnciepdia ID 
map_ID_2_0 = function(x){
  #x = "HSALNG0141642"
  if(!is.character(x)){
    return(print("Error: should enter a string"))
  }
  if(x %in% HSAL_ENSG_lnci[,1]){
    ind_x = which(HSAL_ENSG_lnci[,1] == x)
    return(HSAL_ENSG_lnci[,2][ind_x])
  }else{
    return(NA)
  }
}
#function: HGNC to lncipedia ID 
map_ID_3_0 = function(x){
  if(!is.character(x)){
    return(print("Error: should enter a string"))
  }
  if(x %in% HGNC_ENSG_lnci[,1]){
    ind_x = which(HGNC_ENSG_lnci[,1] == x)
    return(HGNC_ENSG_lnci[,2][ind_x])
  }else{
    return(NA)
  }
}

#function: choose the gene with larger IQR under the same lncip name
#input: count data and lncip_name
#output: lncip_name 
fcn_reduce = function(counts, convert_ID){
  #counts = covid[,-1]
  #lnci_ID = test_2
  #counts = M
  #lnci_ID = lnci_name
  new_ID = convert_ID
  reduc_indx = rep(FALSE, length(new_ID))
  class(counts) = "numeric"
  counts_IQR = apply(counts, 2, iqr)
  for(i in new_ID[!is.na(new_ID)][duplicated(new_ID[!is.na(new_ID)])]){
    #i = new_ID[!is.na(new_ID)][duplicated(new_ID[!is.na(new_ID)])][1]
    comp_ind = which(new_ID == i)
    new_ID[comp_ind[-which.max(counts_IQR[comp_ind])]] = NA
    reduc_indx[comp_ind[-which.max(counts_IQR[comp_ind])]] = TRUE
  }
  return(list(new_ID, reduc_indx))
}


##########################################################################################################################
#function: input a expression matrix for lncRNA with colnames being the original name; output a (reduced) matrix with it's lncipedia name
# input: expression matrix M (row=sample; col=lncRNA)
# type = "ENSG/ENST"/"HSALN"/"HGNC"
# drop: if you want to drop lncRNA with no lncipedia name; otherwise the original names keep remained
# reduce: if you want to reduce the lncRNA to the one with the largest IQR when there is more than one original names corresponds to
# a single same lncipedia name
# output: a new matrix with lncipedia name and index for the dropped genes or reduced genes in the original dataset
name_convert_ncrna = function(M, type, drop = FALSE, reduce = FALSE){
  # M = ncrna
  # type = "HSALN"
  # drop = TRUE
  # reduce = TRUE
  
  ori_name = colnames(M)
  lnci_name = rep(NA, ncol(M))
  if(type == "ENSG/ENST"){
    for(i in 1:length(lnci_name)){
      lnci_name[i] = map_ID_1_0(ori_name[i])
    }
  }
  if(type == "HSALN"){
    for(i in 1:length(lnci_name)){
      lnci_name[i] = map_ID_2_0(ori_name[i])
    }
  }
  if(type == "HGNC"){
    for(i in 1:length(lnci_name)){
      lnci_name[i] = map_ID_3_0(ori_name[i])
    }
  }
  print(sum(is.na(lnci_name)))
  if(drop == FALSE && reduce == FALSE){
    new_M = M
    colnames(new_M) = lnci_name
    colnames(new_M)[is.na(lnci_name)] = ori_name[is.na(lnci_name)]
    na_ind = is.na(lnci_name)
    return(list(new_M, na_ind))
  }
  if(drop == TRUE && reduce == FALSE){
    new_M = M[,!is.na(lnci_name)]
    colnames(new_M) = lnci_name[!is.na(lnci_name)]
    drop_ind = is.na(lnci_name)
    return(list(new_M, drop_ind))
  }
  if(drop == TRUE && reduce == TRUE){
    reduc_name = fcn_reduce(M, lnci_name)
    lnci_name = reduc_name[[1]]
    reduc_ind = reduc_name[[2]]
    new_M = M[,!is.na(lnci_name)]
    colnames(new_M) = lnci_name[!is.na(lnci_name)]
    drop_ind = is.na(lnci_name)
    return(list(new_M, reduc_ind, drop_ind))
  }
  if(drop == FALSE && reduce == TRUE){
    reduc_name = fcn_reduce(M, lnci_name)
    lnci_name = reduc_name[[1]]
    reduc_ind = reduc_name[[2]]
    new_M = M
    colnames(new_M) = lnci_name
    colnames(new_M)[is.na(lnci_name)] = ori_name[is.na(lnci_name)]
    na_ind = is.na(lnci_name)
    return(list(new_M, na_ind, reduc_ind))
  }
}


#######test
#KICH dataset
x <- read.csv("./data/TCGA_KIRP_early_lncRNA_top5per.csv", header=T)
rownames(x) <- x[,1]
x <- x[,-1]
X <- t(x)
KICH = X
test_df = as.matrix(KICH)
colnames(test_df) = KICH[,1]
lnci_KICH = name_convert_ncrna(test_df, "ENSG/ENST", TRUE, FALSE)
dim(lnci_KICH[[1]])
#covid dataset
# d4 <- read.table("SRP253951FeatureCounts_Matirx.tsv",
#                  header=T, sep="\t")
# covid = d4
# test_df = as.matrix(t(covid)[-1,])
# colnames(test_df) = covid[,1]
# lnci_covid = name_convert(test_df, "HSALN", TRUE, TRUE)
# length(lnci_covid[[1]])
# length(grep("HSALN", covid[,1]))
#HGNC dataset











