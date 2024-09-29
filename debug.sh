inputfrag="./output_debug/1-ZLAAO90NB_S7509-S7709.sorted.frag.tsv"
cat ${inputfrag} | grep V350270271L2C002R03901216431 > test.frag.tsv
path_to_fa="/Users/hieunguyen/data/resources/hg19.fa"
inputfrag="test.frag.tsv"
debugdir="debugdir"
mkdir -p ${debugdir};
sampleid="test"

cat $inputfrag |\
    awk '{start=$2 - 1; end= $2 - 1 + 4; name= $5; strand = "+"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    >${debugdir}/${sampleid}.forward_endcoord4bp.bed

cat ${inputfrag} | \
      awk '{start=$3 - 1 - 4; end= $3 - 1; name= $5; strand = "-"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${debugdir}/${sampleid}.reverse_endcoord4bp.bed;

bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${debugdir}/${sampleid}.forward_endcoord4bp.bed | \
      awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${debugdir}/${sampleid}.forward_endmotif4bp.txt
bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${debugdir}/${sampleid}.reverse_endcoord4bp.bed | \
    awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${debugdir}/${sampleid}.reverse_endmotif4bp.txt