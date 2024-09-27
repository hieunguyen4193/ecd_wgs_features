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
path.to.qc.report <- args$input
html.name <- args$batch_name
path.to.save.html.report <- args$output
read_length <- args$read_length

# path.to.qc.report <- "/home/hieunguyen/CRC1382/src/ecd_wgs_features/summary_QC.Rmd"
# html.name <- "test.html"
# path.to.save.html.report <- "/home/hieunguyen/CRC1382/src/ecd_wgs_features/output/QC"
# read_length <- 50



rmarkdown::render(input = path.to.qc.report, 
                  params = list(
                    input = "/home/hieunguyen/CRC1382/src/ecd_wgs_features/output/QC",
                    read_length = read_length
                  ),
                  output_file = html.name, 
                  output_dir = path.to.save.html.report)