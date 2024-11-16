#########################################################################
#########################################################################
### Perform Principal Component Analysis by Multi-Dimensional Scaling ###
#########################################################################
#########################################################################

# load required libraries
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("stringr", quietly = TRUE)) {install.packages("stringr")} 
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("edgeR", quietly = TRUE)) {BiocManager::install("edgeR")} 
if (!require("car", quietly = TRUE)) {install.packages("car")} 

library(stringr) # used for string manipulation
library(edgeR) # package required to generate MDS plot
library(car) # package required for ellipse in MDS plot

# set the working directory to where this script is located
script_path = rstudioapi::getSourceEditorContext()$path
workdir_path = dirname(script_path)
setwd(workdir_path)

# load all data 
COUNTS = read.csv(file.path("..", "..", "Data", "Proteomics", "protein_counts.csv"), header = T)
  rownames(COUNTS) = make.unique(COUNTS$symbol, sep = '.')
  COUNTS$symbol = NULL
  COUNTS[is.na(COUNTS)] = 0
  COUNTS = COUNTS[,c(5:52)]
  COUNTS = COUNTS[,c(1:16)]

metadata = read.csv(file.path("..", "..", "Data", "Proteomics", "metadata.csv"), header = T, row.names = 'X')
  
degs = list()
index = match(colnames(COUNTS),metadata$id)
metadata = metadata[index,]
metadata$samples = colnames(fc_obj)

# identify differentially expressed proteins
y = DGEList(COUNTS, lib.size = colSums(COUNTS),
            norm.factors = calcNormFactors(COUNTS),
            samples = metadata$id,
            group = metadata[['condition']])

# filter lowly expressed genes
keep = filterByExpr(y)
y = y[keep,]
y = calcNormFactors(y)
filteredExpr = cpm(y, log=T)

# format groups for multi-dimensional scaling plot
# set colors for each group
# WT_GF = '#b6d0e2', WT_SPF = 'gray', ASO_GF = '#4682b4', ASO_SPF = 'black'
group = metadata$condition
condition_colors = c(WT_SPF = 'gray',
                     ASO_SPF = 'black',
                     WT_GF = '#b6d0e2',
                     ASO_GF = '#4682b4')
colors = str_replace_all(group, pattern = condition_colors)

# initialize plot parameters
mds = plotMDS(filteredExpr, col = colors, pch = 16)

# get coordinates of all MDS plot points
for(condition in unique(group)) {
  group_points = mds$x[metadata$condition == condition]
  group_points = cbind(group_points, mds$y[metadata$condition == condition])
  center = colMeans(group_points)
  cov_matrix = cov(group_points)
  # set confidence interval ellipse at 80%
  radius = sqrt(qchisq(0.8, df = 2)) 
  
  # get ellipse points
  ellipse_points = car::ellipse(center = center,shape = cov_matrix, radius = radius,
                                type = "l", draw = FALSE)
  car::ellipse(center = center,shape = cov_matrix, radius = radius,
               add = TRUE, col = condition_colors[condition],
               fill = T, fill.alpha = 0.1, center.pch = F)
}

# initialized PDF to save final plot as PDF
output_filepath = file.path("..", "Output", "Proteomics", "MDS_plot_with_ellipses.pdf")
pdf(output_filepath)

# set plot parameters
ylab = paste0('Dimension 2 (',round(mds$var.explained[2]*100),'%)')
xlab = paste0('Dimension 1 (',round(mds$var.explained[1]*100),'%)')
ylim = c(-7.5, 12.5); xlim = c(-12.5, 7.5)

plot(mds$x, mds$y,
     col = colors, pch = 16,
     ylab = ylab, xlab = xlab,
     ylim = ylim, xlim = xlim)
legend("topright",
       legend = unique(group)[c(1,2,3,4)],
       col = unique(colors)[c(1,2,3,4)],
       pch = 16, bty = 'n')

for(condition in unique(group)) {
  group_points = mds$x[metadata$condition == condition]
  group_points = cbind(group_points, mds$y[metadata$condition == condition])
  center = colMeans(group_points)
  cov_matrix = cov(group_points)
  radius = sqrt(qchisq(0.8, df = 2)) # set confidence interval ellipse at 80%
  
  car::ellipse(center = center, shape = cov_matrix, radius = radius,
               add = TRUE, col = condition_colors[condition],
               fill = T, fill.alpha = 0.1, center.pch = F)
}

# save plot as PDF
dev.off()
