inputfrag="/Users/hieunguyen/data/FinaleDB/EE87519.hg19.frag.tsv";
path_to_fa="/Users/hieunguyen/data/resources/hg19_no_chr.fa";
cat ${inputfrag} | \
      awk '{start=$2 - 1; end= $2 - 1 + 4; name= $5; strand = "+"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > test.forward_endcoord4bp.bed;

bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed test.forward_endcoord4bp.bed 