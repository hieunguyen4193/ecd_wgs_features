gc()
rm(list = ls())
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
parser <- ArgumentParser()

parser$add_argument("-s", "--input_short_bam", action="store",
                    help="Path to input preprocessed short BAM files")
parser$add_argument("-l", "--input_long_bam", action="store",
                    help="Path to input preprocessed long BAM files")
parser$add_argument("-f", "--input_full_bam", action="store",
                    help="Path to input preprocessed long BAM files")
parser$add_argument("-o", "--output", action="store",
                    help="Path to save all output")

args <- parser$parse_args()

path.to.short.bam <- args$input_short_bam
path.to.long.bam <- args$input_long_bam
path.to.full.bam <- args$input_full_bam
path.to.save.output <- args$output

sampleid <- str_replace(basename(path.to.full.bam, ".bam", ""))
dir.create(path.to.save.output, showWarnings = FALSE), recursive = TRUE)
#####----------------------------------------------------------------------#####
##### CNA
#####----------------------------------------------------------------------#####
print("Generating CNA features ...")
bin1M <- getBinAnnotations(binSize = 1000, genome = "hg19")
bin1M.grange <- makeGRangesFromDataFrame(
    df = subset(bin1M@data, bin1M@data$use == TRUE), 
    seqnames.field = "chromosome", 
    start.field = "start", 
    end.field = "end", 
    keep.extra.columns = TRUE)

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

output.CNA <- calculate_CNA(bin1M, path.to.full.bam)
write.csv(output.CNA[["CNA"]]@assayData$copynumber %>% as.data.frame(), file.path(path.to.save.output, sprintf("%s.CNA.csv", sampleid)))

#####----------------------------------------------------------------------#####
##### short-to-long ratio in 1M bin
#####----------------------------------------------------------------------#####
print("Generating short-to-long ratio feature")
shortread.counts <- binReadCounts(bin1M, bamfiles = path.to.short.bam[[1]])
shortread.counts <- shortread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
colnames(shortread.counts) <- c("region", "short")

longread.counts <- binReadCounts(bin1M, bamfiles = path.to.long.bam[[1]])
longread.counts <- longread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
colnames(longread.counts) <- c("region", "long")

ratiodf <- merge(shortread.counts, longread.counts, by.x = "region", by.y = "region")
ratiodf$ratio_short_long <- ratiodf$short/ratiodf$long
ratiodf$ratio_short_total <- ratiodf$short/(ratiodf$short + ratiodf$long)

write.csv(ratiodf, file.path(path.to.save.output, sprintf("%s.flen_ratio.csv", sampleid)))
