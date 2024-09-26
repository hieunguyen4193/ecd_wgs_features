
start_time <- Sys.time()

library(QDNAseq)
library(Biobase)
library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(hash)
library(Rsamtools)
library(parallel)
library(GenomicRanges)
library(argparse)
library(testit)

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
# path.to.preprocess.folder <- "/media/hieunguyen/HNSD01/src/ECD_features/input/KABG37"
# path.to.save.output <- file.path("/media/hieunguyen/HNSD01/src/ECD_features/test_output/KABG37")
# gwmd.cores <- 30

parser <- ArgumentParser()

parser$add_argument("-i", "--input", action="store",
                    help="Path to input preprocessed BAM files")
parser$add_argument("-o", "--output", action="store",
                    help="Path to save all output")
# parser$add_argument("-c", "--gwmd_cores", action="store",
#                     help="Number of cores used in calculating GWMD feature")

args <- parser$parse_args()

path.to.preprocess.folder <- args$input
path.to.save.output <- file.path(args$output, sprintf("%s_features", basename(args$input)))
gwmd.cores <- args$gwmd_cores

path.to.save.RDS <- file.path(path.to.save.output, "rds")

#####----------------------------------------------------------------------#####
##### PREPARATION
#####----------------------------------------------------------------------#####
print("Running some preparation steps ...")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
dir.create(path.to.save.RDS, showWarnings = FALSE, recursive = TRUE)

tmpdir <- file.path(path.to.save.output, "tmp")
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

all.paths <- list(
                  bamfile = Sys.glob(file.path(path.to.preprocess.folder, "*.sorted.bam")),
                  em.file = Sys.glob(file.path(path.to.preprocess.folder, "EM", "*.endmotif4bp.txt")),
                  flen.file = Sys.glob(file.path(path.to.preprocess.folder, "flen", "*.flen.txt")),
                  short.bam = Sys.glob(file.path(path.to.preprocess.folder, "short_bam", "*.short.bam")),
                  long.bam = Sys.glob(file.path(path.to.preprocess.folder, "long_bam", "*.long.bam")),
                  short.bam_modify_cutoff = Sys.glob(file.path(path.to.preprocess.folder, "short_bam_modify_cutoff", "*.short.bam")),
                  long.bam_modify_cutoff = Sys.glob(file.path(path.to.preprocess.folder, "long_bam_modify_cutoff", "*.long.bam")),
                  nucleosome = Sys.glob(file.path(path.to.preprocess.folder, "nucleosome", "*.csv")),
                  flen_chr.file = Sys.glob(file.path(path.to.preprocess.folder, "flen_chr", "*.flen.txt")),
                  EM_chr.file = Sys.glob(file.path(path.to.preprocess.folder, "EM_chr", "*EM.txt")),
                  nucleosome_chr.file = Sys.glob(file.path(path.to.preprocess.folder, "nucleosome_chr", "*.csv"))
                  )

for (filetype in names(all.paths)){
  
  if (filetype == "flen_chr.file") {
    print(filetype)
  assert(length(all.paths[[filetype]]) == 22)
  } else if (filetype == "EM_chr.file") {
    print(filetype)
  assert(length(all.paths[[filetype]]) == 22)
  } else {
    print(filetype)
  assert(length(all.paths[[filetype]]) == 1)
  }

}

output <- hash()
output.rds <- hash()

#####----------------------------------------------------------------------#####
##### CNA
#####----------------------------------------------------------------------#####
print("Generating CNA features ...")
bin1M <- getBinAnnotations(binSize = 1000, genome = "hg19")
bin1M.grange <- makeGRangesFromDataFrame(df = subset(bin1M@data, bin1M@data$use == TRUE), seqnames.field = "chromosome", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

calculate_CNA <- function(input.bin, input.path){
  readCounts <- binReadCounts(bins = input.bin, bamfiles = input.path)
  
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="log2")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  function.output <- list(CNA = copyNumbersCalled, readCounts = readCountsFiltered)
  return(function.output)  
}

output.CNA <- calculate_CNA(bin1M, all.paths[["bamfile"]][[1]])

output.rds[["CNA"]] <- output.CNA[["CNA"]]
output.rds[["readCounts"]] <- output.CNA[["readCounts"]]
output[["CNA"]] <- output.CNA[["CNA"]]@assayData$copynumber %>% as.data.frame()

#####----------------------------------------------------------------------#####
##### EM
#####----------------------------------------------------------------------#####
print("Generating EM features ...")
generate_em_feature <- function(input.path){
  emdf <- read.csv(input.path, sep = "\t", header = FALSE, col.names = c("read", "motif")) %>% subset(select = c(motif)) 
  emdf$motif <- toupper(emdf$motif)
  
  emdf.count <- table(emdf$motif) %>% as.data.frame()
  if (nrow(emdf.count) != 0){
    colnames(emdf.count) <- c("motif", "count")
    emdf.count <- subset(emdf.count, grepl("N", emdf.count$motif) == FALSE)
    emdf.count <- emdf.count %>% rowwise() %>%
      mutate(freq = count/sum(emdf.count$count))
    emdf.count <- subset(emdf.count, select = c(motif, freq))
    return(emdf.count)    
  }
}

emdf.count <- generate_em_feature(all.paths[["em.file"]][[1]])
output[["EM"]] <- emdf.count

#####----------------------------------------------------------------------#####
##### FLEN
#####----------------------------------------------------------------------#####
print("Generating FLEN features ...")
generate_flen_feature <- function(input.path){
  flendf <- read.csv(input.path, header = FALSE, col.names = c("flen")) %>%
    rowwise() %>%
    mutate(abs.flen = abs(flen))
  if (nrow(flendf) != 0){
    flen.count <- table(flendf$abs.flen) %>% as.data.frame()
    colnames(flen.count) <- c("size", "count")
    flen.count$size <- as.numeric(as.character(flen.count$size))
    flen.count <- subset(flen.count, flen.count$size >= 50 & flen.count$size <= 350)
    
    flen.count <- flen.count %>% rowwise() %>%
      mutate(freq = count / sum(flen.count$count))
    
    flen.count <- flen.count[order(flen.count$size), ]
    
    output.flendf <- data.frame(size = seq(50, 350))
    output.flendf <- merge(output.flendf, flen.count, by.x = "size", by.y = "size", all.x = TRUE)
    output.flendf[is.na(output.flendf$freq), ] <- 0
    
    output.flendf <- subset(output.flendf, select = c(size, freq, count))
    return(output.flendf)  
  }
}

output.flendf <- generate_flen_feature(all.paths[["flen.file"]][[1]])
output[["flen"]] <- output.flendf

#####----------------------------------------------------------------------#####
##### short-to-long ratio
#####----------------------------------------------------------------------#####
print("Generating short-to-long ratio feature")
shortread.counts <- binReadCounts(bin1M, bamfiles = all.paths[["short.bam"]][[1]])
shortread.counts <- shortread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
colnames(shortread.counts) <- c("region", "short")

longread.counts <- binReadCounts(bin1M, bamfiles = all.paths[["long.bam"]][[1]])
longread.counts <- longread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
colnames(longread.counts) <- c("region", "long")

ratiodf <- merge(shortread.counts, longread.counts, by.x = "region", by.y = "region")
ratiodf$ratio_short_long <- ratiodf$short/ratiodf$long
ratiodf$ratio_short_total <- ratiodf$short/(ratiodf$short + ratiodf$long)

output[["flen_ratio"]] <- ratiodf

#####----------------------------------------------------------------------#####
##### CNA short reads
#####----------------------------------------------------------------------#####
print("Generating short-read CNA features ...")
output.short.CNA <- calculate_CNA(bin1M, all.paths[["short.bam"]][[1]])

output.rds[["CNA_short"]] <- output.short.CNA[["CNA"]]
output[["CNA_short"]] <- output.short.CNA[["CNA"]]@assayData$copynumber %>% as.data.frame()


#####----------------------------------------------------------------------#####
##### NUCLEOSOME
#####----------------------------------------------------------------------#####
# Copy the file
success <- file.copy(all.paths[["nucleosome"]][[1]], file.path(path.to.save.output, "GWfeature_Nucleosome.csv"), overwrite = TRUE)
# Check if the file was copied successfully
if (success) {
  cat("File copied successfully!\n")
} else {
  cat("Failed to copy the file.\n")
}


#####----------------------------------------------------------------------#####
##### NUCLEOSOME FOR EACH REGION
#####----------------------------------------------------------------------#####
# Copy the file
success <- file.copy(all.paths[["nucleosome_chr.file"]][[1]], file.path(path.to.save.output, "GWfeature_Nucleosome_for_each_chr.csv"), overwrite = TRUE)
# Check if the file was copied successfully
if (success) {
  cat("File copied successfully!\n")
} else {
  cat("Failed to copy the file.\n")
}


#####----------------------------------------------------------------------#####
##### short-to-long ratio modify cutoff
#####----------------------------------------------------------------------#####
print("Generating short-to-long ratio feature")
shortread.counts <- binReadCounts(bin1M, bamfiles = all.paths[["short.bam_modify_cutoff"]][[1]])
shortread.counts <- shortread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
colnames(shortread.counts) <- c("region", "short")

longread.counts <- binReadCounts(bin1M, bamfiles = all.paths[["long.bam_modify_cutoff"]][[1]])
longread.counts <- longread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
colnames(longread.counts) <- c("region", "long")

ratiodf <- merge(shortread.counts, longread.counts, by.x = "region", by.y = "region")
ratiodf$ratio_short_long <- ratiodf$short/ratiodf$long
ratiodf$ratio_short_total <- ratiodf$short/(ratiodf$short + ratiodf$long)

output[["flen_ratio_modify_cutoff"]] <- ratiodf


#####----------------------------------------------------------------------#####
##### GENERATE FRAGMENT LENGTH DISTRIBUTION FOR EACH REGION
#####----------------------------------------------------------------------#####
generate_flen_feature <- function(input.path){
  flendf <- read.csv(input.path, header = FALSE, col.names = c("flen")) %>%
    rowwise() %>%
    mutate(abs.flen = abs(flen))
  if (nrow(flendf) != 0){
    flen.count <- table(flendf$abs.flen) %>% as.data.frame()
    colnames(flen.count) <- c("size", "count")
    flen.count$size <- as.numeric(as.character(flen.count$size))
    flen.count <- subset(flen.count, flen.count$size >= 50 & flen.count$size <= 350)
    
    flen.count <- flen.count %>% rowwise() %>%
      mutate(freq = count / sum(flen.count$count))
    
    flen.count <- flen.count[order(flen.count$size), ]
    
    output.flendf <- data.frame(size = seq(50, 350))
    output.flendf <- merge(output.flendf, flen.count, by.x = "size", by.y = "size", all.x = TRUE)
    output.flendf[is.na(output.flendf$freq), ] <- 0
    
    output.flendf <- subset(output.flendf, select = c(size, freq))
    return(output.flendf)  
  }
}


target.flendf <- data.frame(size = seq(50, 350))

for (file_path in all.paths[["flen_chr.file"]]) {
  # Extract the region from the file path using a regular expression
  region <- sub(".+\\.(chr\\.[^.]+)\\..+", "\\1", file_path)
  
  region.flendf <- generate_flen_feature(file_path)
  
  if (is.null(region.flendf) == FALSE){
    colnames(region.flendf) <- c("size", region)
    target.flendf <- merge(target.flendf, region.flendf, by.x = "size", by.y = "size", all.x = TRUE)      
  }
}

output[["flen_chr"]] <- target.flendf


#####----------------------------------------------------------------------#####
##### GENERATE EM DISTRIBUTION FOR EACH REGION
#####----------------------------------------------------------------------#####
all.motifs <- c()
for (i in c("A", "T", "G", "C")){
  for (j in c("A", "T", "G", "C")){
    for (k in c("A", "T", "G", "C")){
      for (l in c("A", "T", "G", "C")){
        all.motifs <- c(all.motifs, sprintf("%s%s%s%s", i, j, k, l))
      }
    }
  }
}

target.emdf <- data.frame(motif = all.motifs)

for (file_path in all.paths[["EM_chr.file"]]) {
  # Extract the region from the file path using a regular expression
  region <- sub(".+\\.(chr\\.[^.]+)\\..+", "\\1", file_path)
  
  region.emdf <- generate_em_feature(file_path)
  
  if (is.null(region.emdf) == FALSE){
    colnames(region.emdf) <- c("motif", region)
    target.emdf <- merge(target.emdf, region.emdf, by.x = "motif", by.y = "motif", all.x = TRUE)      
  }
}

output[["EM_chr"]] <- target.emdf


##### CLEAN UP

print("Saving features to output dir")
for (feat in names(output)){
  write.csv(output[[feat]], file.path(path.to.save.output, sprintf("GWfeature_%s.csv", feat)))
}

for (feat in names(output.rds)){
  saveRDS(output.rds[[feat]], file.path(path.to.save.output, sprintf("GWfeature_%s.rds", feat)))
}

end_time <- Sys.time()
start_time - end_time
