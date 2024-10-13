#!/bin/bash
#SBATCH --job-name="download_finaleddb"
#SBATCH --partition=research_hi
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --output=/mnt/NAS_PROJECT/vol_ECDteam/hieunho/data/finaledb/log.log


# cd /mnt/NAS_PROJECT/vol_ECDteam/hieunho/data/finaledb

# for i in {87786..88325};
# do
    # wget https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/entries/EE${i}/hg19/EE${i}.hg19.frag.tsv.bgz
    # wget https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/entries/EE${i}/hg19/EE${i}.hg19.insert_size_metrics.txt
# done



cd /mnt/NAS_PROJECT/vol_ECDteam/hieunho/data/finaledb

# Path to the CSV file
csv_file=/mnt/NAS_PROJECT/vol_ECDteam/hieunho/data/finaledb/metadata_frag_tsv.csv

# Loop through the CSV file, skipping the header
tail -n +2 "$csv_file" | while IFS=',' read -r SampleID frag_tsv_file
do
    name=${frag_tsv_file##*/}
    if [ -e "$name" ]; then
        echo "Skip" $name
    else
        wget https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/${frag_tsv_file}
    fi
done

