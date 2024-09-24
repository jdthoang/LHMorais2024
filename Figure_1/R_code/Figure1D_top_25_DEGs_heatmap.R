###################################################################
###################################################################
### Generate Scaled Heatmap of Top 25 DEGs for Each Comparison ###
###################################################################
###################################################################

# load required packages
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")} 
if (!require("ComplexHeatmap", quietly = TRUE)) {BiocManager::install("ComplexHeatmap")}
if (!require("viridis", quietly = TRUE)) {install.packages("viridis")}

library("DESeq2")
library("ComplexHeatmap")
library("viridis")

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

CPM = sweep(COUNTS,2,as.numeric(colSums(COUNTS)),FUN="/")*1000000

# diffential gene expression with DESeq2
dds = DESeq(DESeqDataSetFromMatrix(countData=COUNTS, 
                             colData=metadata, 
                             design=~condition, tidy = F))

# set conditions to compare
comp = c('ASO_GF','WT_GF')
comp = rbind(comp,c('ASO_SPF', 'WT_SPF'))
comp = rbind(comp,c('WT_GF', 'WT_SPF'))
comp = rbind(comp,c('ASO_GF', 'ASO_SPF'))
             
# create heatmaps of top 25 DEG
for(i in 1:4){
  con1 = comp[i,1]; con2 = comp[i,2]; logfcthr = 0.5
  res = results(dds, contrast = c('condition',con1, con2), pAdjustMethod ='BH'); res = res[order(res$padj),]
  topDEGs = res[(!is.na(res$padj)),] #res$padj < 0.05 & 
  topDEGs = topDEGs[order(topDEGs$log2FoldChange, decreasing = T),]
  
  # select genes to plot
  sig_genes = c(head(rownames(topDEGs), 25),tail(rownames(topDEGs), 25))
  
  # convert to matrix and scale
  mat = CPM[sig_genes, metadata$condition %in% c(con1, con2)]
  heat = t(scale(t(mat)))
  
  # sort by metadata
  metadata$Genotype = toupper(metadata$Genotype)
  ann = data.frame(
    Genotype = metadata$Genotype[metadata$condition %in% c(con1, con2)],
    Microbiome = metadata$Microbiome[metadata$condition %in% c(con1, con2)],
    stringsAsFactors = FALSE)
  
  # create the colour mapping
  colours = list(Genotype = c('WT' = 'white', 'ASO' = 'black'), Microbiome = c('SPF' = 'gray', 'GF' = '#4682b4'))
  
  # create heatmap structure
  colAnn = HeatmapAnnotation(
    df = ann,
    which = 'col',
    na_col = 'white',
    col = colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    border = T,
    annotation_legend_param = list(
      Genotype = list(
        nrow = 2, title = 'Genotype', title_position = 'topcenter', legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold')),
      Microbiome = list(
        nrow = 2, title = 'Microbiome', title_position = 'topcenter', legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold'))))
  
  hmap = Heatmap(heat, cluster_rows = T, cluster_columns = F,
                 top_annotation = colAnn, width = unit(0.3, "npc"), height = unit(0.75, "npc"),
                 heatmap_legend_param = list(at = seq(-4, 4, 2)),
                 rect_gp = gpar(col = "black"),
                 col = viridis(100))
  
  # initialize SVG to save file
  filename = paste0(con1, '_vs_', con2, '_top-25-DEG-heatmap.svg',sep = "")
  filepath = file.path("..", "Output", "Transcriptomics", filename)
  svg(file=filepath,
      width = 3.47*3, height = 6.27*2)
  draw(hmap, heatmap_legend_side = 'left', annotation_legend_side = 'right', row_sub_title_side = 'left',
       padding = unit(c(0.5, 0.5, 0.5, 0.5), "in"))
  dev.off()
}
