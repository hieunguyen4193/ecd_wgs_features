  # Sort the reference BED file
  nucleosome_ref=/Users/hieunguyen/data/resources/rpr_map_EXP0779.bed
  sort -k 1V,1 -k 2n,2 ${nucleosome_ref} -o ${nucleosome_ref%.bed*}.sorted.bed