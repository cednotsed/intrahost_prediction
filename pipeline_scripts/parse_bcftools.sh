wkdir=/mnt/c/git_repos/intrahost_prediction
ref=$wkdir/data/genomes/MN908947.3.fna
in_dir=$wkdir/results/pipeline_out.020225/vcf_out
out_dir=$wkdir/results/pipeline_out.020225/vcf_out.biallelic

in_dir=$wkdir/results/pipeline_out.tonkin/vcf_out
out_dir=$wkdir/results/pipeline_out.tonkin/vcf_out.biallelic
#biosample=SAMD00268090

biosample=$1

in_vcf=$in_dir/$biosample.bcftools.vcf.gz
out_vcf=$out_dir/$biosample.bcftools.biallelic.vcf.gz

bcftools norm \
    -f $ref \
    -m -any \
    -o $out_vcf \
    -O z \
    $in_vcf

