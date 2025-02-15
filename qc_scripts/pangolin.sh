genome_path=../data/alignments/reassembled.gap_filtered.aln
genome_path=../data/alignments/tonkin.gap_filtered.aln
out_path=$(echo $genome_path|sed "s|.aln|.pango_lineage.csv|g"|sed "s|alignments|metadata|g")

pangolin \
    --threads 12 \
    --outfile $out_path \
    $genome_path


