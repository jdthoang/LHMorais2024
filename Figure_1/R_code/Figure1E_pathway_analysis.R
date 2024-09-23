#########################################################
#########################################################
### Identification of Striatum Enriched Gene Pathways ###
#########################################################
#########################################################

# load required libraries
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")}
if (!require("RITAN", quietly = TRUE)) {BiocManager::install("RITAN")}
if (!require("RITANdata", quietly = TRUE)) {BiocManager::install("RITANdata")}
if (!require("tidyverse", quietly = TRUE)) {install.packages("tidyverse")}
if (!require("ggplot2", quietly = TRUE)) {install.packages("ggplot2")}
if (!require("ggtree", quietly = TRUE)) {BiocManager::install("ggtree")}

library(DESeq2)
library(RITAN)
library(RITANdata)
library(tidyverse)
library(ggplot2)
library(ggtree)

# set WD
path = "path"; setwd(path)

# read in count matrix and run DESEQ
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

# setup comparisons to run
comp = c('ASO_SPF','WT_SPF')
comp = rbind(comp,c('ASO_GF', 'WT_GF'))
comp = rbind(comp,c('WT_GF', 'WT_SPF'))
comp = rbind(comp,c('ASO_GF', 'ASO_SPF'))

# set pathway databases to use for analysis
resources = c('KEGG_filtered_canonical_pathways','GO_slim_generic','ReactomePathways')

up_terms = data.frame()
down_terms = data.frame()
logfcthr = 0.5; pthr = 0.05
# run pathway analysis to initialize number of unique terms
for(i in 1:4){
  # find list of DEGs
  con1 = comp[i,1]; con2 = comp[i,2]
  res = results(dds, contrast = c('condition',con1, con2), pAdjustMethod ='BH'); res = res[order(res$padj),]
  DEGs = res[(res$padj < pthr & !is.na(res$padj)),]
  # find upregulated and downregulated terms
  upreg = rownames(DEGs[DEGs$log2FoldChange > logfcthr,])
  downreg = rownames(DEGs[DEGs$log2FoldChange < -logfcthr,])
  e_upreg = term_enrichment(toupper(upreg), resources = resources); e_upreg = e_upreg[e_upreg$q < pthr,]
    e_upreg$GeneRatio = e_upreg$n/e_upreg$n.set; e_upreg = e_upreg[order(e_upreg$GeneRatio, decreasing = T),]
  e_downreg = term_enrichment(toupper(downreg), resources = resources); e_downreg = e_downreg[e_downreg$q < pthr,]
    e_downreg$GeneRatio = e_downreg$n/e_downreg$n.set; e_downreg = e_downreg[order(e_downreg$GeneRatio, decreasing = T),]
  up_terms = rbind(up_terms,head(e_upreg,10))
  down_terms = rbind(down_terms,head(e_downreg,10))
  }

# create dataframes to input values
GOs = unique(c(up_terms$name, down_terms$name)); nGOs = length(GOs); print(nGOs)
data_up = data.frame("GOs" = rep(GOs, 4), 
                   "Condition" = rep(paste0(comp[,1],' over ',comp[,2]), each = nGOs),
                   "GeneRatio" = '', 
                   "p.adjust" = '')
data_down = data.frame("GOs" = rep(GOs, 4), 
                     "Condition" = rep(paste0(comp[,1],' over ',comp[,2]), each = nGOs),
                     "GeneRatio" = '', 
                     "p.adjust" = '')

# rerun and input data into dataframe
for(i in 1:4){
  con1 = comp[i,1]; con2 = comp[i,2]
  res = results(dds, contrast = c('condition',con1, con2), pAdjustMethod ='BH'); res = res[order(res$padj),]
  DEGs = res[(res$padj < 0.05 & !is.na(res$padj)),]
  upreg = rownames(DEGs[DEGs$log2FoldChange > logfcthr,])
  downreg = rownames(DEGs[DEGs$log2FoldChange < -logfcthr,])
  e_upreg = term_enrichment(toupper(upreg), resources = resources)
    e_upreg$GeneRatio = e_upreg$n/e_upreg$n.set
      e_upreg = e_upreg[order(e_upreg$GeneRatio, decreasing = T),]
  e_downreg = term_enrichment(toupper(downreg), resources = resources)
    e_downreg$GeneRatio = e_downreg$n/e_downreg$n.set
      e_downreg = e_downreg[order(e_downreg$GeneRatio, decreasing = T),]
  
  index = match(data_up$GOs[(1:nGOs)+nGOs*(i-1)],e_upreg$name)
  data_up$GeneRatio[(1:nGOs)+nGOs*(i-1)] = e_upreg$n[index]/e_upreg$n.set[index]
  data_up$p.adjust[(1:nGOs)+nGOs*(i-1)] = as.numeric(e_upreg$q[index])

  index = match(data_down$GOs[(1:nGOs)+nGOs*(i-1)],e_downreg$name)
  data_down$GeneRatio[(1:nGOs)+nGOs*(i-1)] = e_downreg$n[index]/e_downreg$n.set[index]
  data_down$p.adjust[(1:nGOs)+nGOs*(i-1)] = as.numeric(e_downreg$q[index])}

# clean up data
data_up$GeneRatio = as.numeric(data_up$GeneRatio)
data_up$p.adjust = as.numeric(data_up$p.adjust)
data_up$GeneRatio[data_up$GeneRatio == 0] = NA
data_down$GeneRatio = as.numeric(data_down$GeneRatio)
data_down$p.adjust = as.numeric(data_down$p.adjust)
data_down$GeneRatio[data_down$GeneRatio == 0] = NA

# log transform q values
data_up$logq = -log10(data_up$p.adjust)
data_down$logq = -log10(data_down$p.adjust)

# make data square to calculate euclidean distance
mat_up = data_up %>% 
  select(-logq, -p.adjust) %>%
  pivot_wider(names_from = Condition, values_from = GeneRatio) %>% 
  data.frame()
mat_up[is.na(mat_up)] = 0
row.names(mat_up) = mat_up$GOs
mat_up = mat_up[,-1]

mat_down = data_down %>% 
  select(-logq, -p.adjust) %>%
  pivot_wider(names_from = Condition, values_from = GeneRatio) %>% 
  data.frame()
mat_down[is.na(mat_down)] = 0
row.names(mat_down) = mat_down$GOs
mat_down = mat_down[,-1]

mat = cbind(mat_up, mat_down)
clust = hclust(dist(mat %>% as.matrix()))

# create dendrogram
ddgram = as.dendrogram(clust)
ggtree_plot = ggtree(ddgram)

data_up$Condition = paste0("UP     ",data_up$Condition)
data_down$Condition = paste0("DOWN     ",data_down$Condition)
data = rbind(data_up, data_down)

data$GOs = factor(data$GOs, levels = clust$labels[clust$order])
data$Condition = factor(data$Condition, levels = unique(data$Condition)[c(1,2,3,4,5,6,7,8)])

data$GOs = gsub('GO_slim_generic.','GO: ',data$GOs)
data$GOs = gsub('ReactomePathways.','Reactome: ',data$GOs)
data$GOs = gsub('KEGG_filtered_canonical_pathways.KEGG','KEGG: ',data$GOs)
data$GOs = gsub('_',' ',data$GOs)
data$GOs = gsub('/.','',data$GOs)

dotplot = ggplot(data = data, aes(x=Condition, y = GOs, color = logq, size = GeneRatio)) + 
  geom_vline(aes(xintercept = Condition), color = "lightgray", alpha = 0.5) +
  geom_hline(aes(yintercept = GOs), color = "lightgray", alpha = 0.5) +
  geom_point() + 
  scale_color_viridis_c(name = '-log10(FDR)') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_y_discrete(position = "right")
