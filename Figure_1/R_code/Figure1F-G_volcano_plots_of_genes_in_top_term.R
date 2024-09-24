################################################
################################################
### Volcano Plot of Genes in Top Term (Mito) ###
################################################
################################################

# load required packages
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")}
if (!require("RITANdata", quietly = TRUE)) {BiocManager::install("RITANdata")}
if (!require("EnhancedVolcano", quietly = TRUE)) {BiocManager::install("EnhancedVolcano")}

library("DESeq2") # used for differential gene expression analysis
library("RITANdata") 
library("EnhancedVolcano")

# set the working directory to where this script is located
script_path = rstudioapi::getSourceEditorContext()$path
workdir_path = dirname(script_path)
setwd(workdir_path)

# read in data
COUNTS = read.table(file.path("..", "..", "Data", "Transcriptomics", "gene_counts.tsv"), 
                    sep = '\t', 
                    header = T,
                    row.names = 'X')
                    
metadata = read.csv(file.path("..", "..", "Data", "Transcriptomics", "metadata.csv"), 
                    row.names = 'X')

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
  filename = paste0(con1,"_vs_",con2,"_volcano.png")
  filepath = file.path("..", "Output", "Transcriptomics", filename)
  ggsave(filepath, units = 'in', width = 6, height = 8)
}
