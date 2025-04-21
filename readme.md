# Order of scripts

## Metadata parsing
### Parse BioSample metadata
parse_ncbi_biosample_xml_to_csv.full.ipynb
parse_biosample_meta.early.R
parse_biosample_meta.transitions.R
get_biosample_submission_date.sh

### Parse GISAID metadata
filter_gisaid_metadata.R
merge_early_and_transition_metadata.R

## Bioinfo pipeline
### Process samples
consensus_pipeline.V3.transition.sbs

### Generate consensus
generate_consensus_from_vcf.V2.R
parse_bcftools.sh # Multiallelic->biallelic VCFs

### Remove gappy sequences
merge_and_filter_alignment.R

### Pangolin lineage assignment
pangolin.sh

### Remove misassigned genomes and SNP count outliers
calculate_snp_distances.R
pango_filter.R

## GISAID data
### Calculate monthly SAV frequencies
submit_monthly_freq_script.sbs
calculate_monthly_frequencies.geo.crick.R # Split by region
calculate_monthly_frequencies.V2.crick.R  # All time

### Aggregate SAV fitness
aggregate_monthly_prop.R
aggregate_monthly_prop.all_time.R
aggregate_monthly_prop.cross_dataset.R

## Interhost linkage
### Generate presence absence matrix of mutations for D' calculations
split_gisaid_by_month.R
generate_presence_matrix.V2.crick.R
submit_generate.crick.sbs

### Calculate D'
calculate_Dprime.parallel.V2.crick.R
submit_Dprime.crick.sbs



