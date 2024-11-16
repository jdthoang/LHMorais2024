###############################################
###############################################
### Volcano Plot of Proteins in Top 5 Terms ###
###############################################
###############################################

# load required packages
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("RITANdata", quietly = TRUE)) {BiocManager::install("RITANdata")}
if (!require("EnhancedVolcano", quietly = TRUE)) {BiocManager::install("EnhancedVolcano")}

library(RITANdata)
library(EnhancedVolcano)

# set the working directory to where this script is located
script_path = rstudioapi::getSourceEditorContext()$path
workdir_path = dirname(script_path)
setwd(workdir_path)

# load mitocarta data
url <- 'https://personal.broadinstitute.org/scalvo/MitoCarta_Download/MitoPathways3.0.gmx'
mitopathways = read.delim(file = url)
geneset_list$MitoCarta = mitopathways

# find and read in files with differentially expressed proteins
filepath = file.path("..", "..", "Data", "Proteomics")
files = list.files(path = filepath, 
                   pattern = 'dep')

logfcthr = 0.5; pthr = 0.05
for(i in 1:4){
  dep = read.csv(file.path(filepath, (files[i])))
  dep$gene_symbol = make.unique(dep$gene_symbol,'.')
  
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
  both_mito = both_mito[c('Detoxification','ROS_and_glutathione_metabolism')]
  both_mito = unique(unlist(both_mito)); both_mito = both_mito[both_mito != ""]
  both_mito = both[toupper(both) %in% both_mito]
  
  # Create a filtered subset of the data for labeling (including Snca and both_mito)
  label_genes = c(both_mito, "Timm10", "Timm9", "Tomm22", "Dars2", "Tars2", "Vars2", "Sars2", "Nars2")
  selected_genes = dep[dep$gene_symbol %in% label_genes, , drop = FALSE]
  
  # Extract the colors from keyvals to match selected_genes
  selected_genes$label_color = keyvals[match(selected_genes$gene_symbol, dep$gene_symbol)]
  
  # Split the selected_genes by color category
  blue_genes = selected_genes[selected_genes$label_color == "blue",]
  red_genes = selected_genes[selected_genes$label_color == "red",]
  
  dep$p_value
  # Base volcano plot without the default labels
  base_plot = EnhancedVolcano(dep,
                              title = paste0(con2, '(-) vs ', con1, '(+)'), subtitle = '',
                              lab = dep$gene_symbol, selectLab = c(''),
                              x = 'log2_foldchange', FCcutoff = logfcthr, xlab = bquote(bold(~log[2] ~ 'fold change')), xlim = c(-2.5, 2.5),
                              y = 'p_value', ylab = bquote(bold(~-log[10] ~ bold(P))), pCutoff = pthr, ylim = c(0, 10),
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
                     aes(x = log2_foldchange, y = -log10(p_value), label = gene_symbol),
                     fill = "white", color = "blue", size = 8, fontface = 'bold',
                     box.padding = 0.5, point.padding = 0.5,
                     segment.color = 'black', min.segment.length = 0.1, segment.size = 0.5,
                     direction = "both", xlim = c(NA, -0.5), ylim = c(1, NA), 
                     max.time = 60, max.iter = 100000,
                     force = 100, force_pull = 100,
                     max.overlaps = 50)
  
  # Red labels
  base_plot = base_plot +
    geom_label_repel(data = red_genes,
                     aes(x = log2_foldchange, y = -log10(p_value), label = gene_symbol),
                     fill = "white", color = "red", size = 8, fontface = 'bold',
                     box.padding = 0.5, point.padding = 0.5,
                     segment.color = 'black', min.segment.length = 0.1, segment.size = 0.5,
                     direction = "both", xlim = c(0.5, NA), ylim = c(1, NA),
                     max.time = 60, max.iter = 100000,
                     force = 100, force_pull = 100,
                     max.overlaps = 50)
  
  # Darkgrey labels

  
  # Print the final plot
  print(base_plot)
  
  # save volcano plot
  filetype = '.png'
  filename = paste0(con1,"_vs_",con2,"_volcano_rm nonsig",filetype)
  imgpath = file.path("..", "Output", "Proteomics", filename)
  ggsave(imgpath, units = 'in', width = 13, height = 10)
}
