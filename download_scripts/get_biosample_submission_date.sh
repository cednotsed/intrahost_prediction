biosample_list=../data/metadata/all_sra_metadata.biosamples_only.txt
file_path=../data/metadata/all_sra_metadata.submission_dates.tsv

export NCBI_API_KEY='ed474d19fd845a6bc6f1b0dba7e803a17308'

echo 'BioSample\tSubmissionDate' \
> $file_path

cat $biosample_list | \
    xargs -n 1 \
    sh -c 'echo "$0";
           efetch -db biosample -id "$0" -format xml| \
               xtract -pattern BioSampleSet \
               -element BioSample@accession BioSample@submission_date \
           >> '"$file_path"''

