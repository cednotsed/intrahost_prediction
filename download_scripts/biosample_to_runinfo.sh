acc_list=../data/external_datasets/tonkin/tonkin.biosamples.txt
out_path=../data/external_datasets/tonkin/runinfo.tonkin.csv

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
