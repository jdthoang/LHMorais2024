################################################
################################################
### Volcano Plot of Genes in Top Term (Mito) ###
################################################
################################################

# load required packages
library(DESeq2)
library(RITANdata)
library(EnhancedVolcano)

# set WD
path = "path"; setwd(path)

# read in data
COUNTS = read.table("./gene_counts_striatum.tsv", sep = '\t', header = T)
metadata = read.csv("./metadata_striatum.csv", row.names = 'X')

# diffential gene expression with DESeq2
dds = DESeq(DESeqDataSetFromMatrix(countData=COUNTS, 
                             colData=metadata, 
                             design=~condition, tidy = F))

# set conditions to compare
comp = c('ASO_GF','WT_GF')
comp = rbind(comp,c('ASO_SPF', 'WT_SPF'))
comp = rbind(comp,c('WT_GF', 'WT_SPF'))
comp = rbind(comp,c('ASO_GF', 'ASO_SPF'))
             
logfcthr = 0.5; pthr = 0.05
# Create the contrast matrix for DE
for(i in 1:4){
  con1 = comp[i,1]; con2 = comp[i,2]
  res = results(dds, contrast = c('condition',con1, con2), pAdjustMethod ='BH'); res = res[order(res$padj),]
  
  # set colors for each group
  keyvals = ifelse((res$log2FoldChange < -logfcthr & res$padj < pthr), 'blue',
                   ifelse((res$log2FoldChange > logfcthr & res$padj < pthr), 'red',
                          'darkgrey')); keyvals[is.na(keyvals)] = 'darkgrey'
  names(keyvals)[keyvals == 'red'] = con1
  names(keyvals)[keyvals == 'darkgrey'] = 'unchanged'
  names(keyvals)[keyvals == 'blue'] = con2
  
  # subset genes that are up and downregulated
  up = c(rownames(res)[res$padj < pthr & res$log2FoldChange > logfcthr]); up = up[!is.na(up)]
  down = c(rownames(res)[res$padj < pthr & res$log2FoldChange < -logfcthr]); down = down[!is.na(down)]
  both = c(up,down); print(paste0(con1,' vs ',con2,' - ',length(both),' DEGs'))
  
  # filter those genes for mitochondrion involved genes according to GO terms
  both_mito = unlist(geneset_list$ReactomePathways$`Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.`)
  both_mito = both_mito[both_mito != ""]
  both_mito = both[toupper(both) %in% both_mito]
  
  # make volcano plot
  p = EnhancedVolcano(res,
                  title = paste0(con2,'(-) vs ',con1,'(+)'), subtitle = '',
                  lab = rownames(res), selectLab = c(both_mito,'Snca'),
                  x = 'log2FoldChange', FCcutoff = logfcthr, xlim = c(-8,8),
                  y = 'padj', ylab = bquote(~-Log[10] ~ italic(FDR)), pCutoff = pthr, ylim = c(0,10),
                  pointSize = 1.5,
                  labSize = 7.0, labCol = 'black', labFace = 'bold', boxedLabels = T,
                  colAlpha = 0.8, colCustom = keyvals,
                  legendPosition = 'top', legendLabSize = 14, legendIconSize = 4.0,
                  drawConnectors = TRUE, widthConnectors = 0.5, typeConnectors = 'open',  colConnectors = 'black',
                  max.overlaps = 100, min.segment.length = 0.1, arrowheads = F,
                  cutoffLineCol = 'black')
  
  # save volcano plot
  ggsave(paste0(con1,"_vs_",con2,"_volcano.png"), units = 'in', width = 6, height = 8)
}