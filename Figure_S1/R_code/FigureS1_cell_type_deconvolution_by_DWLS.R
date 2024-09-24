#################################################################################
#################################################################################
### Cell Type Composition Deconvolution Using Dampened Weighted Least Squares ###
#################################################################################
#################################################################################

# load required packages
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("rstudioapi", quietly = TRUE)) {install.packages("rstudioapi")} 
if (!require("MAST", quietly = TRUE)) {BiocManager::install("MAST")}
if (!require("DWLS", quietly = TRUE)) {install.packages("DWLS")} 
if (!require("Seurat", quietly = TRUE)) {install.packages("Seurat")} 

library(DWLS)
library(Seurat)

# load sc data and metadata
sc_path = 'path to scRNAseq matrix files'
data = ReadMtx(mtx = paste0(scpath, '/matrix.mtx'),
               cells = paste0(scpath, './barcodes.tsv'),
               features = paste0(scpath, './genes.tsv'))
data = as.matrix(data)
metadata = read.csv('metadata.csv')

# create cell type signatures
label = as.vector(metadata$cell_type)
signature = buildSignatureMatrixMAST(scdata = data,
                                     id = label,
                                     path = './DWLS.results/',
                                     diff.cutoff = 0.5,
                                     pval.cutoff = 0.01)
save(signature,file="DWLS.RData")

signature = MASTSignatureMatrixGivenDE(scdata = data,
                                       id = labels,
                                       path = 'DWLS.results',
                                       diff.cutoff = 0.5,
                                       pval.cutoff = 0.01)
save(signature,file="./DWLS.results/DWLS.RData")

# read in count matrix for deconvolution
COUNTS = read.csv("count_matrix.tsv", header = T)
CPM = sweep(COUNTS,2,as.numeric(colSums(COUNTS)),FUN="/")*1000000

# filter counts for genes in signature matrix
CPM.trim = CPM[rownames(signature),]
samples = colnames(CPM.trim)

# Initializing Results Data Frames for Each Algorithm that DWLS employs
DWLS = data.frame()

# Computing Cell Fractions for Each Sample and Adding to Results Data Frames
for(sample in samples){
  bulk = CPM.trim[,sample]
  names(bulk) = rownames(CPM.trim)
  tr=trimData(signature,bulk)
  solDWLS=solveDampenedWLS(tr$sig,tr$bulk)
  DWLS = rbind(DWLS,solDWLS)}
model = list(DWLS,SVR,OLS)

# Adding the rownames and columns back to each dataframe in list
L = lapply(model, function(df){colnames(df) = colnames(signature); rownames(df) = samples; df})
write.csv(L[1],file="CellFractions_PFCref_DWLS.csv")