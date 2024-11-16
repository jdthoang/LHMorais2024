################################################
################################################
### Volcano Plot of Genes in Top Term (Mito) ###
################################################
################################################

# Load required packages
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")}
if (!require("RITANdata", quietly = TRUE)) {BiocManager::install("RITANdata")}
if (!require("EnhancedVolcano", quietly = TRUE)) {BiocManager::install("EnhancedVolcano")}

library("DESeq2") # Used for differential gene expression analysis
library("RITANdata") 
library("EnhancedVolcano")

# Set the working directory to where this script is located
script_path = rstudioapi::getSourceEditorContext()$path
workdir_path = dirname(script_path)
setwd(workdir_path)

# Read in data
COUNTS = read.table(file.path("..", "..", "Data", "Transcriptomics", "gene_counts.tsv"), 
                    sep = '\t', 
                    header = T)
  rownames(COUNTS) = make.unique(COUNTS$genes); COUNTS$genes = NULL                    
metadata = read.csv(file.path("..", "..", "Data", "Transcriptomics", "metadata.csv"), 
                    row.names = 'X')

# Diffential gene expression with DESeq2
dds = DESeq(DESeqDataSetFromMatrix(countData=COUNTS,
                             colData=metadata, 
                             design=~condition, tidy = F))

# Set conditions to compare
comp = c('ASO_GF','WT_GF')
comp = rbind(comp,c('ASO_SPF', 'WT_SPF'))
comp = rbind(comp,c('WT_GF', 'WT_SPF'))
comp = rbind(comp,c('ASO_GF', 'ASO_SPF'))
             
logfcthr = 0.5; pthr = 0.05
# Create the contrast matrix for DE
for(i in 1:4){
  con1 = comp[i,1]; con2 = comp[i,2]
  res = results(dds, contrast = c('condition',con1, con2), alpha = 0.05, pAdjustMethod = 'BH')
  res = res[order(res$padj),]
  
  # Set colors for each group
  keyvals = ifelse((res$log2FoldChange < -logfcthr & res$padj < pthr), 'blue',
                   ifelse((res$log2FoldChange > logfcthr & res$padj < pthr), 'red',
                          'darkgrey')); keyvals[is.na(keyvals)] = 'darkgrey'
  names(keyvals)[keyvals == 'red'] = con1
  names(keyvals)[keyvals == 'darkgrey'] = 'unchanged'
  names(keyvals)[keyvals == 'blue'] = con2
  
  # Subset genes that are up and downregulated
  up = c(rownames(res)[res$padj < pthr & res$log2FoldChange > logfcthr]); up = up[!is.na(up)]
  down = c(rownames(res)[res$padj < pthr & res$log2FoldChange < -logfcthr]); down = down[!is.na(down)]
  both = c(up,down); print(paste0(con1,' vs ',con2,' - ',length(both),' DEGs'))
  
  # Filter those genes for mitochondrion involved genes according to GO terms
  both_mito = unlist(geneset_list$ReactomePathways$`Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.`)
  both_mito = both_mito[both_mito != ""]
  both_mito = both[toupper(both) %in% both_mito]
  
  # Create a filtered subset of the data for labeling (including Snca and both_mito)
  label_genes = c(both_mito, 'Snca' ,'hSNCA')
  selected_genes = res[rownames(res) %in% label_genes, , drop = FALSE]
  
  # Extract the colors from keyvals to match selected_genes
  selected_genes$label_color = keyvals[match(rownames(selected_genes), rownames(res))]
  
  # Split the selected_genes by color category
  blue_genes = selected_genes[selected_genes$label_color == "blue", ]
  darkgrey_genes = selected_genes[selected_genes$label_color == "darkgrey", ]
  red_genes = selected_genes[selected_genes$label_color == "red", ]
  
  # Base volcano plot without the default labels
  base_plot = EnhancedVolcano(res,
                              title = paste0(con2, '(-) vs ', con1, '(+)'), subtitle = '',
                              lab = rownames(res), selectLab = c(''),
                              x = 'log2FoldChange', FCcutoff = logfcthr, xlab = bquote(bold(~log[2] ~ 'fold change')), xlim = c(-8, 8),
                              y = 'padj', ylab = bquote(bold(~-log[10] ~ bold(FDR))), pCutoff = pthr, ylim = c(0, 10),
                              axisLabSize = 25, borderWidth = 1.5,
                              pointSize = 1.5,
                              labSize = 7.0, boxedLabels = FALSE, parseLabels = FALSE,
                              colAlpha = 0.8, colCustom = keyvals,
                              legendPosition = 'top', legendLabSize = 14, legendIconSize = 4.0,
                              drawConnectors = FALSE,
                              max.overlaps = 100, min.segment.length = 0.1,
                              cutoffLineCol = 'black')
  
  # Add labels manually for each color group with appropriate limits
  # Blue labels
  base_plot = base_plot +
    geom_label_repel(data = blue_genes,
                     aes(x = log2FoldChange, y = -log10(padj), label = rownames(blue_genes)),
                     fill = "white", color = "blue", size = 8, fontface = 'bold',
                     box.padding = 0.5, point.padding = 0.5,
                     segment.color = 'black', min.segment.length = 0.2, segment.size = 0.5,
                     direction = "both", xlim = c(NA, -0.5), ylim = c(1, NA), 
                     max.time = 60, max.iter = 100000,
                     force = 100, force_pull = 50,
                     max.overlaps = 50)
  
  # Red labels
  base_plot = base_plot +
    geom_label_repel(data = red_genes,
                     aes(x = log2FoldChange, y = -log10(padj), label = rownames(red_genes)),
                     fill = "white", color = "red", size = 8, fontface = 'bold',
                     box.padding = 0.5, point.padding = 0.5,
                     segment.color = 'black', min.segment.length = 0.2, segment.size = 0.5,
                     direction = "both", xlim = c(0.5, NA), ylim = c(1, NA),
                     max.time = 60, max.iter = 100000,
                     force = 100, force_pull = 50,
                     max.overlaps = 50)
  
  # Darkgrey labels
  base_plot = base_plot +
    geom_label_repel(data = darkgrey_genes,
                     aes(x = log2FoldChange, y = -log10(padj), label = rownames(darkgrey_genes)),
                     fill = "white", color = "darkgrey", size = 6, fontface = 'bold',
                     box.padding = 1, point.padding = 1,
                     segment.color = 'black',min.segment.length = 0.2, segment.size = 0.5,
                     direction = "both", ylim = c(0, 0.5),
                     force = 200, force_pull = 0.05,
                     max.overlaps = 100)
  
  # Print the final plot
  print(base_plot)
  
  # save volcano plot
  filetype = '.png'
  filename = paste0(con1,"_vs_",con2,"_volcano",filetype)
  filepath = file.path("..", "Output", "Transcriptomics", filename)
  ggsave(filepath, units = 'in', width = 13, height = 10)
}
