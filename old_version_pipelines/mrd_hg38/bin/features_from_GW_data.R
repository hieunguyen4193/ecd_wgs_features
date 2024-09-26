#!/usr/bin/env Rscript
# coding: utf-8

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

parser <- ArgumentParser()

parser$add_argument("-i", "--input", action="store",
                    help="Path to input preprocessed BAM files")
parser$add_argument("-s", "--sample", action="store",
                    help="Labcode")
parser$add_argument("-o", "--output", action="store",
                    help="Path to save all output")
parser$add_argument("-f", "--feature_type", action="store",
                    help="Type of feature")
                    

args <- parser$parse_args()

sample <- args$sample
feature_type <- args$feature_type

path.to.preprocess.folder <- args$input
path.to.save.output <- args$output


output <- hash()

#####----------------------------------------------------------------------#####
##### EM
#####----------------------------------------------------------------------#####
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

                colnames(emdf.count)[which(colnames(emdf.count) == "freq")] <- sample
                return(emdf.count)
              }
}

#####----------------------------------------------------------------------#####
##### CNA
#####----------------------------------------------------------------------#####
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

#####----------------------------------------------------------------------#####
##### FLEN
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
                colnames(output.flendf)[which(colnames(output.flendf) == "freq")] <- sample
                return(output.flendf)  
              }
}

#####----------------------------------------------------------------------#####
##### PREPARATION
#####----------------------------------------------------------------------#####
print("Running some preparation steps ...")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

if (feature_type == 'EM') {

            all.paths <- list(EM = Sys.glob(file.path(path.to.preprocess.folder, sprintf("%s.EM.txt", sample))))

            for (filetype in names(all.paths)){
              print(filetype)
              assert(length(all.paths[[filetype]]) == 1)
            }

            print("Generating EM features ...")
            emdf.count <- generate_em_feature(all.paths[["EM"]][[1]])
            output[["EM"]] <- emdf.count

} else if (feature_type == 'CNA_SHORT') {

            all.paths <- list(CNA_SHORT = Sys.glob(file.path(path.to.preprocess.folder, sprintf("%s.short.bam", sample))))

            for (filetype in names(all.paths)){
              print(filetype)
              assert(length(all.paths[[filetype]]) == 1)
            }

            print("Generating CNA_SHORT features ...")
            bin1M <- getBinAnnotations(binSize = 1000, genome = "hg19")
            bin1M.grange <- makeGRangesFromDataFrame(df = subset(bin1M@data, bin1M@data$use == TRUE), seqnames.field = "chromosome", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

            output.short.CNA <- calculate_CNA(bin1M, all.paths[["CNA_SHORT"]][[1]])

            output[["CNA_SHORT"]] <- output.short.CNA[["CNA"]]@assayData$copynumber %>% as.data.frame()

} else if (feature_type == 'FLEN_RATIO') {

            all.paths <- list(short.bam = Sys.glob(file.path(path.to.preprocess.folder, sprintf("%s.short.bam", sample))),
                              long.bam = Sys.glob(file.path(path.to.preprocess.folder, sprintf("%s.long.bam", sample))))

            for (filetype in names(all.paths)){
              print(filetype)
              assert(length(all.paths[[filetype]]) == 1)
            }

            print("Generating short-to-long ratio feature")
            bin1M <- getBinAnnotations(binSize = 1000, genome = "hg19")
            bin1M.grange <- makeGRangesFromDataFrame(df = subset(bin1M@data, bin1M@data$use == TRUE), seqnames.field = "chromosome", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

            shortread.counts <- binReadCounts(bin1M, bamfiles = all.paths[["short.bam"]][[1]])
            shortread.counts <- shortread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
            colnames(shortread.counts) <- c("region", "short")

            longread.counts <- binReadCounts(bin1M, bamfiles = all.paths[["long.bam"]][[1]])
            longread.counts <- longread.counts@assayData$counts %>% as.data.frame() %>% rownames_to_column("region")
            colnames(longread.counts) <- c("region", "long")

            ratiodf <- merge(shortread.counts, longread.counts, by.x = "region", by.y = "region")
            ratiodf$ratio_short_long <- ratiodf$short/ratiodf$long
            ratiodf$ratio_short_total <- ratiodf$short/(ratiodf$short + ratiodf$long)
            colnames(ratiodf)[which(colnames(ratiodf) == "ratio_short_long")] <- sprintf("%s.ratio_short_long", sample)
            colnames(ratiodf)[which(colnames(ratiodf) == "ratio_short_total")] <- sprintf("%s.ratio_short_total", sample)

            output[["FLEN_RATIO"]] <- ratiodf

} else if (feature_type == 'CNA') {

            all.paths <- list(bamfile = Sys.glob(file.path(path.to.preprocess.folder, sprintf("%s.sorted.bam.sorted.bam", sample))))

            for (filetype in names(all.paths)){
              print(filetype)
              assert(length(all.paths[[filetype]]) == 1)
            }

            print("Generating CNA features ...")
            bin1M <- getBinAnnotations(binSize = 1000, genome = "hg19")
            bin1M.grange <- makeGRangesFromDataFrame(df = subset(bin1M@data, bin1M@data$use == TRUE), seqnames.field = "chromosome", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

            output.CNA <- calculate_CNA(bin1M, all.paths[["bamfile"]][[1]])

            output[["CNA"]] <- output.CNA[["CNA"]]@assayData$copynumber %>% as.data.frame()

} else if (feature_type == 'FLEN') {

            all.paths <- list(FLEN = Sys.glob(file.path(path.to.preprocess.folder, sprintf("%s.FLEN.txt", sample))))

            for (filetype in names(all.paths)){
              print(filetype)
              assert(length(all.paths[[filetype]]) == 1)
            }

            print("Generating FLEN features ...")
            output.flendf <- generate_flen_feature(all.paths[["FLEN"]][[1]])
            output[["FLEN"]] <- output.flendf

}


##### CLEAN UP

print("Saving features to output dir")
for (feat in names(output)){
  write.csv(output[[feat]], file.path(path.to.save.output, sprintf("%s_GWfeature_%s.csv", sample, feat)))
}

end_time <- Sys.time()
start_time - end_time
