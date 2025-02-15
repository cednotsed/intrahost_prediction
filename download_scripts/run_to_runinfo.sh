#acc_list=../data/metadata/sra_metadata/filtered_sra_accessions.early.txt
#out_path=/mnt/c/git_repos/wuhu_rooting/data/metadata/sra_metadata/runinfo.early.csv

acc_list=../data/metadata/sra_metadata/filtered_sra_accessions.early.missing.txt
out_path=../data/metadata/sra_metadata/runinfo.early.missing.csv

acc_list=../data/metadata/sra_metadata/filtered_sra_accessions.transition.missing.txt
out_path=../data/metadata/sra_metadata/runinfo.transition.missing.csv

# Get headers
esearch \
    -db sra \
    -query ERS3350306 \
    | efetch \
    -format runinfo \
    |head -n1 \
> $out_path

while read acc <&3
do
    echo $acc

    esearch \
        -db sra \
        -query $acc \
        | efetch \
        -format runinfo \
        |tail -n +2 \
        >> $out_path
done 3< $acc_list
