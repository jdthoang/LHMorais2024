###############################################
###############################################
### Volcano Plot of Proteins in Top 5 Terms ###
###############################################
###############################################

# load required packages
if (!require("RITANdata", quietly = TRUE)) {BiocManager::install("RITANdata")}
if (!require("EnhancedVolcano", quietly = TRUE)) {BiocManager::install("EnhancedVolcano")}

library(RITANdata)
library(EnhancedVolcano)

# set WD
path = "path"; setwd(path)

# load mitocarta data
url <- 'https://personal.broadinstitute.org/scalvo/MitoCarta_Download/MitoPathways3.0.gmx'
mitopathways = read.delim(file = url)
geneset_list$MitoCarta = mitopathways

# find and read in files with differentially expressed proteins
filepath = file.path("..", "..", "Data", "Proteomics")
files = list.files(path = filepath, 
                   pattern = 'dep')

logfcthr = 0.5; pthr = 0.05
# Create the contrast matrix for DE


for(i in 1:4){
  dep = read.csv(file.path(filepath, (files[i])))
  
  con1 = unlist(strsplit(files[i], '[_.]+'))[2]
  con2 = unlist(strsplit(files[i], '[_.]+'))[4]
  
  # set colors for each group
  keyvals = ifelse((dep$log2_foldchange < -logfcthr & dep$p_value < pthr), 'blue',
                   ifelse((dep$log2_foldchange > logfcthr & dep$p_value < pthr), 'red',
                          'darkgrey')); keyvals[is.na(keyvals)] = 'darkgrey'
  names(keyvals)[keyvals == 'red'] = con1
  names(keyvals)[keyvals == 'darkgrey'] = 'unchanged'
  names(keyvals)[keyvals == 'blue'] = con2
  
  # subset genes that are up and downregulated
  up_dep = c(dep$gene_symbol[dep$p_value < pthr & dep$log2_foldchange > logfcthr]); up_dep = up_dep[!is.na(up_dep)]
  down_dep = c(dep$gene_symbol[dep$p_value < pthr & dep$log2_foldchange < -logfcthr]); down_dep = down_dep[!is.na(down_dep)]
  both = c(up_dep,down_dep)
  
  # filter those genes for mitochondrion involved genes according to GO terms
  both_mito = geneset_list$MitoCarta
  both_mito = both_mito[c('Detoxification','ROS_and_glutathione_metabolism',
                          'mt.tRNA_synthetases', 'Protein_import_and_sorting',
                          'Lipid_metabolism')]
  both_mito = unique(unlist(both_mito)); both_mito = both_mito[both_mito != ""]
  both_mito = both[toupper(both) %in% both_mito]
  
  # make volcano plot
  EnhancedVolcano(dep,
                  title = paste0(con2,'(-) vs ',con1,'(+)'), subtitle = '', subtitleLabSize = 0,
                  lab = dep$gene_symbol, selectLab = both_mito,
                  x = 'log2_foldchange', FCcutoff = 0.5, xlim = c(-3,3),
                  y = 'p_value', ylab = bquote(~-Log[10] ~ italic(P)), pCutoff = 5e-2, ylim = c(0,10),
                  pointSize = 1.5,
                  labSize = 7.0, labCol = 'black', labFace = 'bold', boxedLabels = T,
                  colAlpha = 0.8, colCustom = keyvals,
                  legendPosition = 'top', legendLabSize = 14, legendIconSize = 4.0,
                  drawConnectors = TRUE, widthConnectors = 0.5, typeConnectors = 'open',  colConnectors = 'black',
                  max.overlaps = 100, min.segment.length = 0.1, arrowheads = F,
                  cutoffLineCol = 'black',
                  axisLabSize = 27)
  
  # save volcano plot
  filename = file.path("..", "Output", "Proteomics", paste0('protein_',con1,"_vs_",con2,".pdf"))
  ggsave(filename, units = 'in', width = 12.5, height = 10.5)
}
