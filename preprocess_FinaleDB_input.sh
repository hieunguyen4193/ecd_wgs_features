while getopts "i:o:f:r:n:c:" opt; do
  case ${opt} in
    i )
      inputfrag=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] input .frag file [-o] outputdir"
      exit 1
      ;;
  esac
done

sampleid=$(echo ${inputfrag} | xargs -n 1 basename)
sampleid=${sampleid%.tsv*}
mkdir -p ${outputdir};

echo -e "Working on the file " ${sampleid} " at the link " ${inputfrag}

# make the following modification on the column of FinaleDB input file
# 1. Add "chr" to the first column
# 2. Add a new column at $5 (readID) by concatenating the first 5 columns.
# 3. Move column 4 (read mappping quality) to column 6
# 4. Replace column 4 (read length) with the fragment length (column 3 - column 2)
# 5. Filter out rows where the second column (fragment start) is less than 5, avoid NaN in the output in calculating End motif 4bp
# 6. Keep only rows where the first column is in 1-22, X, Y, or MT

cat ${inputfrag} | awk -v OFS='\t' '{$1="chr"$1;$5 = $1"_"$2"_"$3"_"$4"_"$5;$6 = $4;$4 = $3 - $2; if ($2 > 5){print $0}}' > ${outputdir}/${sampleid}.tmp.frag.tsv
awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/' ${outputdir}/${sampleid}.tmp.frag.tsv > ${outputdir}/${sampleid}.FinaleDB.frag.tsv
rm -rf ${inputfrag%.frag*}.tmp.frag.tsv