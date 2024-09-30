gc()
rm(list = ls())

library(argparse)

parser <- ArgumentParser()

parser$add_argument("-i", "--input", action="store",
                    help="Path to input QC folder containing all samples in the batch")
parser$add_argument("-o", "--output", action="store",
                    help="Path to save all output")
parser$add_argument("-n", "--html_name", action="store",
                    help="Name of the output html file")
parser$add_argument("-r", "--read_length", action="store",
                    help="Path to save output html file")

args <- parser$parse_args()

path.to.qc.report <- args$input # path to the input QC folder.
html.name <- args$batch_name # name of the output html file.
path.to.save.html.report <- args$output # path to save the output html file.
read_length <- args$read_length # read length of the sequencing data. Default: set 50.

read_length <- as.integer(read_length)
##### path to the Rmd file
path.to.rmd.file <- "/Users/hieunguyen/src/ecd_wgs_features/summary_QC.Rmd"

rmarkdown::render(input = path.to.rmd.file, 
                  params = list(
                    input = path.to.qc.report,
                    read_length = read_length
                  ),
                  output_file = html.name, 
                  output_dir = path.to.save.html.report)