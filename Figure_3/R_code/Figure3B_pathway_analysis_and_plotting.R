############################################################
############################################################
### Identification of Striatum Enriched Protein Pathways ###
############################################################
############################################################

# load required libraries
library(RITAN)
library(RITANdata)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(cowplot)

# find files with differentially expressed proteins
filepath = file.path("..", "..", "Data", "Proteomics")
files = list.files(path = filepath, 
                   pattern = '.csv')

# read in mitocarta pathways and load into RITAN
mitopathways = read.delim("path to mitocarta")
geneset_list$MitoCarta = mitopathways
resources = c('MitoCarta')

terms = data.frame()
logfcthr = 0.5; pthr = 0.05
# run pathway analysis to initialize number of unique terms
for(i in 1:4){
  # find list of DEPs
  res = read.csv(files[i]); res = res[order(res$p_value),] 
  res$log2_foldchange = res$log2_foldchange*flip[i]
  DEPs = res[(res$p_value < pthr & !is.na(res$p_value)),]
  DEPs = DEPs[!is.na(DEPs$gene_symbol),]
  rownames(DEPs) = make.unique(DEPs$gene_symbol)
  
  # find DEPs with logfcthr cutoff
  DEPs = rownames(DEPs[abs(DEPs$log2_foldchange) > logfcthr,])
  e = term_enrichment(toupper(DEPs), resources = resources)
    e = e[e$p < pthr,]
    e$GeneRatio = e$n/e$n.set; e = e[order(e$GeneRatio, decreasing = T),]
  terms = rbind(terms,head(e,10))
  }

terms$name[duplicated(terms$name)]

# create dataframes to input values
GOs = unique(terms$name); nGOs = length(GOs); print(nGOs)
data = data.frame("GOs" = rep(GOs, 4), 
                  "Condition" = rep(files, each = nGOs),
                  "GeneRatio" = '', 
                  "pval" = '')

# rerun and input data into dataframe
for(i in 1:4){
  # find list of DEPs
  res = read.csv(files[i]); res = res[order(res$p_value),] 
  res$log2_foldchange = res$log2_foldchange*flip[i]
  DEPs = res[(res$p_value < pthr & !is.na(res$p_value)),]
  DEPs = DEPs[!is.na(DEPs$gene_symbol),]
  rownames(DEPs) = make.unique(DEPs$gene_symbol)
  
  # find enriched terms
  DEPs = rownames(DEPs[abs(DEPs$log2_foldchange) > logfcthr,])
  e = term_enrichment(toupper(DEPs), resources = resources)
    e = e[e$p < pthr,]
    e$GeneRatio = e$n/e$n.set; e = e[order(e$GeneRatio, decreasing = T),]
  
  index = match(data$GOs[(1:nGOs)+nGOs*(i-1)],e$name)
  data$GeneRatio[(1:nGOs)+nGOs*(i-1)] = e$n[index]/e$n.set[index]
  data$pval[(1:nGOs)+nGOs*(i-1)] = as.numeric(e$p[index])
}

# clean up data
data$GeneRatio = as.numeric(data$GeneRatio)
  data$pval = as.numeric(data$pval)
  data$GeneRatio[data$GeneRatio == 0] = NA
  data$GOs = gsub('MitoCarta.','',data$GOs)
  data$GOs = gsub('_',' ',data$GOs)
  
# log transform p values
data$logp = -log10(data$pval)

# make data square to calculate euclidean distance
mat = data %>% 
  select(-GeneRatio, -pval) %>%
  pivot_wider(names_from = Condition, values_from = logp) %>% 
  data.frame()
mat[is.na(mat)] = 0
row.names(mat) = mat$GOs
mat = mat[,-1]

clust = hclust(dist(mat %>% as.matrix()))

# create dendrogram
ddgram = as.dendrogram(clust)
ggtree_plot = ggtree(ddgram)

data$GOs = factor(data$GOs, levels = clust$labels[clust$order])
data$Condition = factor(data$Condition, levels = unique(data$Condition)[c(1,2,3,4)])


output_filepath = file.path("..", "Output", "Proteomics", "proteomics_pathways_both.pdf")
pdf(output_filepath, width = 3.47*2.3, height = 6.27*2)
dotplot = ggplot(data = data, aes(x=Condition, y = GOs, color = logp, size = GeneRatio)) + 
  geom_vline(aes(xintercept = Condition), color = "lightgray", alpha = 0.5) +
  geom_hline(aes(yintercept = GOs), color = "lightgray", alpha = 0.5) +
  geom_point() + 
  scale_color_viridis_c(name = '-log10(P)') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_y_discrete(position = "right")

plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.2,-0.05, 2), align = 'h')
dev.off()
