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

cat ${inputfrag} | awk -v OFS='\t' '{$5 = $1"_"$2"_"$3"_"$4"_"$5;$6 = $4;$4 = $3 - $2; if ($2 > 5){print $0}}' > ${outputdir}/${sampleid}.tmp.frag.tsv

# Filter rows where the first column is in 1-22, X, Y, or MT
awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/' ${outputdir}/${sampleid}.tmp.frag.tsv > ${outputdir}/${sampleid}.FinaleDB.frag.tsv
rm -rf ${inputfrag%.frag*}.tmp.frag.tsv