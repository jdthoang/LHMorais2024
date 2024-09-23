#######################################################################
#######################################################################
### Align Raw Data to Custom Alpha-Syn Overexpressing Mouse Genome ###
#######################################################################
#######################################################################

# load required packages
library(Rsubread)

# set Working directory
path = "path"; setwd(path)

# get list of fastq files from working directory and name list of outputs
fastq.files = c(list.files(path = './', pattern = ".fastq.gz$",
                         full.names = TRUE, recursive = TRUE)); fastq.files
bam.files = paste0(substr(basename(fastq.files),1,nchar(basename(fastq.files))-9),".BAM")
 
# build the index
Rsubread::buildindex(basename="aso_mm_39", reference="./aso_mm.fa.gz", memory = 45000)

# align the reads to reference and generate BAM files
Rsubread::align(index = "aso_mm_39",
                readfile1 = fastq.files,
                type = "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                nthreads = 18,
                output_file = bam.files,
                unique = F, nBestLocations = 1)

# convert aligned reads to feature level counts
bam.files = c(list.files(path = './', pattern = ".BAM$",
                           full.names = TRUE, recursive = TRUE)); bam.files
fc = Rsubread::featureCounts(files = bam.files, annot.inbuilt = NULL, annot.ext = './custom_gtf/aso_mm.gtf',
                   GTF.featureType = 'exon',
                   GTF.attrType = 'gene_id',
                   GTF.attrType.extra = 'gene_name',
                   isGTFAnnotationFile = TRUE,
                   isPairedEnd = FALSE, nthreads = 18)

# save gene counts file as tsv
COUNTS = fc$counts
  rownames(COUNTS) = make.unique(COUNTS$row.names, sep = '.')
  COUNTS$row.names = NULL
write.table(COUNTS, "./gene_counts_striatum.tsv", sep = '\t')