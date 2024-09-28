bash 01_BAM_to_FRAGS.sh \
    -i /Users/hieunguyen/data/tmp/debug_ecd_wgs_features/data/1-ZLAAO90NB_S7509-S7709.sorted.bam \
    -o ./output_debug/ \
    -n 10 \
    -f /Users/hieunguyen/data/resources/hg19.fa;

bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh \
    -i ./output_debug/1-ZLAAO90NB_S7509-S7709.sorted.frag.tsv  \
    -o ./output_debug/ \
    -f /Users/hieunguyen/data/resources/hg19.fa \
    -r /Users/hieunguyen/data/resources/rpr_map_EXP0779.sorted.bed;

