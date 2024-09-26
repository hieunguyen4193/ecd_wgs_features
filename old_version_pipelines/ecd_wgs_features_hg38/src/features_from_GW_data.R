
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
                  em.file = Sys.glob(file.path(path.to.preprocess.folder, "EM", "*.endmotif4bp.txt")),
                  flen.file = Sys.glob(file.path(path.to.preprocess.folder, "flen", "*.flen.txt")),
                  nucleosome = Sys.glob(file.path(path.to.preprocess.folder, "nucleosome", "*.csv"))
                  )

for (filetype in names(all.paths)){
  print(filetype)
  assert(length(all.paths[[filetype]]) == 1)
}

output <- hash()
output.rds <- hash()

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
