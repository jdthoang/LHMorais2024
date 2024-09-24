############################################################
############################################################
### Identification of Striatum Enriched Protein Pathways ###
############################################################
############################################################

# load required libraries
if (!require("RITAN", quietly = TRUE)) {BiocManager::install("RITAN")} 
if (!require("RITANdata", quietly = TRUE)) {BiocManager::install("RITANdata")} 
if (!require("tidyverse", quietly = TRUE)) {install.packages("tidyverse")} 
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("ggplot2", quietly = TRUE)) {install.packages("ggplot2")} 
if (!require("ggtree", quietly = TRUE)) {BiocManager::install("ggtree")} 
if (!require("cowplot", quietly = TRUE)) {install.packages("cowplot")} 

library(RITAN); library(RITANdata)
library(tidyverse)
library(ggplot2); library(ggtree)
library(cowplot)

# set the workind directory to the main directory of the repository
script_path = rstudioapi::getSourceEditorContext()$path
workdir_path = dirname(script_path)
setwd(workdir_path)


# find files with differentially expressed proteins
filepath = file.path("..", "..", "Data", "Proteomics")
files = list.files(path = filepath, 
                   pattern = 'dep')

# read in mitocarta pathways and load into RITAN
url <- 'https://personal.broadinstitute.org/scalvo/MitoCarta_Download/MitoPathways3.0.gmx'
mitopathways = read.delim(file = url)
geneset_list$MitoCarta = mitopathways
resources = c('MitoCarta')

terms = data.frame()
logfcthr = 0.5; pthr = 0.05
# run pathway analysis to initialize number of unique terms
for(i in 1:4){
  # find list of DEPs
  res = read.csv(file.path(filepath, (files[i]))); res = res[order(res$p_value),] 
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
  res = read.csv(file.path(filepath, (files[i]))); res = res[order(res$p_value),] 
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
