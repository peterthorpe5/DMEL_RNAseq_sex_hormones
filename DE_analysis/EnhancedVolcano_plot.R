#Volcano plot with EnhancedVolcano (https://github.com/kevinblighe/EnhancedVolcano)
#  BiocManager::install('EnhancedVolcano')
#with results<- a dataframe containing the gene ids, log2FoldChange, adjusted-pvalues (padj)
#label_volcano<-a dataframe with the gene ids to label.




setwd("C:/Users/pjt6/OneDrive - University of St Andrews/Bat_FOXP2_overexpression_Dec_2021_final/volcano_plot")
library(EnhancedVolcano)
library(magrittr)
in_data = read.table("data.txt", header=T,stringsAsFactors=F)
genes = in_data$Gene
padj = in_data$padj
# simple  log2FoldChange	padj
 EnhancedVolcano(in_data,
    lab = rownames(in_data),
    x = 'log2FoldChange',
    y = 'pvalue', 
    xlim = c(-5,5), 
    ylim = c(0,4),
    title = 'GFP_striatum_vs_KI_striatum', 
    pCutoff = 0.05, 
    FCcutoff = 1.2, 
    pointSize = 1.5, labSize = 2.0,
    legendLabSize = 10,
    colAlpha = 1)
    
  
# genes as labels  
 EnhancedVolcano(in_data,
    lab = genes,
    x = 'log2FoldChange',
    y = 'pvalue', 
    xlim = c(-5,5), 
    ylim = c(0,4),
    title = 'GFP_striatum_vs_KI_striatum', 
    pCutoff = 0.05, 
    FCcutoff = 1.2, 
    pointSize = 1.5, labSize = 2.0,
    legendLabSize = 10,
    colAlpha = 1)
    

 EnhancedVolcano(in_data,
    lab = "",
    x = 'log2FoldChange',
    y = 'pvalue', 
    xlim = c(-5,5), 
    ylim = c(0,4),
    title = 'GFP_striatum_vs_KI_striatum', 
    pCutoff = 0.05, 
    FCcutoff = 1.2, 
    pointSize = 1.5, labSize = 2.0,
    legendLabSize = 10,
    colAlpha = 1)

   # example 
# EnhancedVolcano(in_data, lab=in_data$gene, selectLab = unique((label_volcano$Gene)), x = 'log2FoldChange', y = 'padj',xlim = c(-2.5, 2.5),xlab = bquote(~Log[2]~ 'fold change'),ylab = bquote(~-Log[10]~adjusted~italic(P)),FCcutoff = 0.5,pCutoff = 0.001,colAlpha = 0.6,pointSize = c(ifelse(results$Gene%in%Sec$Gene, 3, 1)),gridlines.major = FALSE,gridlines.minor = FALSE,vline=0.5,boxedLabels = FALSE,labCol = "black",labFace = "bold",title = "Regulated genes",shapeCustom = keyvals.shape_sec)
