gc()
rm(list = ls())

path.to.qc.report <- "/home/hieunguyen/CRC1382/src/ecd_wgs_features/summary_QC.Rmd"
html.name <- "test.html"
path.to.save.html.report <- "/home/hieunguyen/CRC1382/src/ecd_wgs_features/output/QC"
read_length <- 50
rmarkdown::render(input = path.to.qc.report, 
                  params = list(
                    input = "/home/hieunguyen/CRC1382/src/ecd_wgs_features/output/QC",
                    read_length = read_length
                  ),
                  output_file = html.name, 
                  output_dir = path.to.save.html.report)