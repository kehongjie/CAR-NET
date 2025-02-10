# test------
#gene_list <- read.csv("pathway analysis/TCGA_KIRP_early_gene_top5per.csv",header=T)
# x <- c("TP53", "EGFR", "MYC", "CDKN2A", "PTEN")
# background <- c("TP53", "EGFR", "MYC", "CDKN2A", "PTEN", "BRCA1", "BRCA2", "KRAS", "NRAS", "BRAF")

load("./data/pathways.rda")
library(ggplot2)


## Main function to make pathway enrichment plot
run.path <- function(gene.list, background, path.interest=NULL) {
  ## select pathways of interest
  if(is.null(path.interest)) {
    path.interest <- c(GOBP.genesets, GOCC.genesets, GOMF.genesets, KEGG.genesets, 
                       Reactome.genesets, Biocarta.genesets)
  }
  
  path_run <- path.fisher(x=gene.list, background=background,
                          pathway=path.interest)
  if (nrow(path_run$pval)>=1) {
    path_result <- path_run$pval
    path_result$pathway <- gsub(",", ";", path_result$pathway)
    path_result$gene_count <- sapply(strsplit(path_result$gene, ";"), length)  
    path_result$GeneRatio <- path_result$gene_count / path_result$pathway_size
    path_result$p_adjust <- path_result$qvalue  
    path_result$log_p_adjust <- -log10(path_result$p_adjust)
    
    color_gradient <- scale_color_gradient(low = "#dedddd", high = "#e3768e")
    
    plt <- ggplot(path_result, aes(x = GeneRatio, y = pathway, size = gene_count, color = log_p_adjust)) +
      geom_point(alpha = 0.8) +  
      theme_bw(base_size=15) +
      scale_size(range = c(3, 10)) +  
      color_gradient +  
      labs(x = "Gene ratio", y = "Pathway name",
           size = "Count", color = "-log10(p-value)") +
      theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "right",
        panel.grid.major = element_line(color = "gray80", linetype = "dashed")
      )
    print(plt)
  } else {
    print("NO significant pathway!")
  }

}


######### internal functions for run.path ########## 

path.fisher <- function(x, background, pathway, d=0) {
  if(d==0) {d <- length(pathway)}
  
  x<-toupper(x)
  background<-toupper(background)
  pathway <- lapply(pathway,toupper)
  pathway <- lapply(pathway,function(a) intersect(background,a))
  back_non_sel <- setdiff(background, x)
  
  
  result <- matrix(NA, length(pathway), 5)
  match_gene_list <- vector("list", length(pathway))
  colnames(result) <- c("pvalue", "Selected_in_path", "Selected_not_in_path",
                        "NonSele_in_path", "NonSele_not_in_path")
  rownames(result) <- names(match_gene_list) <- names(pathway)
  
  for (i in 1:length(pathway)) {
    path <- unlist(pathway[i])
    count_table <- matrix(0,2,2)
    
    ####in the gene list and in the pathway
    count_table[1,1] <- length(intersect(x, path))
    
    ####in the gene list but not in the pathway
    count_table[1,2] <- length(x)-count_table[1,1]
    
    ####not in the gene list but in the pathway
    count_table[2,1] <- length(intersect(back_non_sel, path))
    
    ####not in the gene list and not in the pathway
    count_table[2,2] <- length(back_non_sel)-count_table[2,1]       
    
    pval <- fisher.test(count_table, alternative="greater")$p
    
    match_gene_list[[i]] <- intersect(x, path)
    # match_num<-length(matched_gene)
    
    result[i,] <- c(pval, as.vector(t(count_table)))
  }
  
  order_by_pval <- order(result[,"pvalue"], decreasing=FALSE)
  result <- result[order_by_pval, ]
  match_gene_list <- match_gene_list[order_by_pval]
  
  path_size <- as.numeric(result[,"Selected_in_path"]+result[,"NonSele_in_path"])
  pos_size <- which(path_size>=5 & path_size<=100)
  pos_pval <- which(result[,"pvalue"]<0.05)
  # pos_pval <- order(result[,"pvalue"])[1:10] ## temporary 
  result <- result[intersect(pos_size,pos_pval),]
  match_gene_list <- match_gene_list[intersect(pos_size,pos_pval)]
  result[,"pvalue"] <- round(result[,"pvalue"], 5)
  # result <- result[1:d,]
  # match_gene_list <- match_gene_list[1:d]
  
  qvalue <- round(p.adjust(result[,"pvalue"], "BH"), 5)
  # result <- cbind(result, qvalue)
  result <- data.frame(pathway=rownames(result), result, 
                       pathway_size=path_size[intersect(pos_size,pos_pval)], 
                       qvalue=qvalue, 
                       gene=as.vector(sapply(match_gene_list, function(x) paste(x,collapse=";"))))
  rownames(result) <- NULL
  
  # pval_table <- cbind(rownames(result), result)
  # colnames(pval_table)[1] <- "pathway"
  
  # return(list(pval=result[1:d,], match_gene=match_gene_list[1:d]))
  return(list(pval=result, match_gene=match_gene_list))
}
