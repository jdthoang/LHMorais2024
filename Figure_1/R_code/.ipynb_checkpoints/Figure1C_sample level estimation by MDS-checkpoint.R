#########################################################################
#########################################################################
### Perform Principal Component Analysis by Multi-Dimensional Scaling ###
#########################################################################
#########################################################################

# load required libraries
library(DESeq2) # used for differential gene expression analysis
library(stringr) # used for string manipulation
library(edgeR) # package required to generate MDS plot
library(car) # package required for ellipse in MDS plot

# set WD
path = "path"; setwd(path)

# read in count matrix generated from alignment
COUNTS = read.table("./gene_counts.tsv", sep = '\t', header = T)
metadata = read.csv("./metadata.csv", row.names = 'X')

# differential gene expression with DESeq2
dds = DESeq(DESeqDataSetFromMatrix(countData = COUNTS,
                             colData = metadata,
                             design = ~condition, tidy = F))

# format groups for multi-dimensional scaling plot
# set colors for each group
# WT_GF = '#b6d0e2', WT_SPF = 'gray', ASO_GF = '#4682b4', ASO_SPF = 'black'
group = as.character(dds$condition)
condition_colors = c(WT_SPF = 'gray',
                     ASO_SPF = 'black',
                     WT_GF = '#b6d0e2',
                     ASO_GF = '#4682b4')
colors = str_replace_all(group, pattern = condition_colors)

# initialize plot parameters
mds = plotMDS(dds, col = colors, pch = 16)

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
pdf(paste0(path,'MDS plot with ellipses.pdf'))

# set plot parameters
ylab = paste0(mds$axislabel,' 2 (',round(mds$var.explained[2]*100),'%)')
xlab = paste0(mds$axislabel,' 1 (',round(mds$var.explained[1]*100),'%)')
ylim = c(-1.5, 2.5); xlim = c(-2.5, 2.5)

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