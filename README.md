# de novo tools

Scripts dedicated to de novo short variants (SNVs & indels)

## DeNovoMiner.py
This script runs simple bcftools filtration steps to identify de novo variant candidates from a multi-sample vcf including child-father-mother trios.
It has been tested on Deepvariant 1.5 / Glnexus variants calls obtained on Illumina WGS, and an example of config file for this case study is provided.
<br>
Usage: python DeNovoMiner.py -v <multi_sample_vcf_file> -p <pedigree_file> -c <config_yaml> -o <output_dir>
<br>
Manual reviewing of candidates can be performed on the variant calls using :
https://github.com/francois-lecoquierre/genomics_shortcuts/blob/main/classify_vcf_from_igv_v1.0.py

## DNM_cluster_by_sample.py
Takes as input a multisample vcf file of de novo mutations and identifies clusters.
<br>
Usage: python DNM_cluster_by_sample.py -v <vcf_file> -o <output_folder> -c <max_distance_between_two_var_to_call_cluster>

## DNM_distance_between_samples.py
Takes as input a multisample vcf file of de novo mutations and identifies variants at close proximity between samples
<br>
Usage: python DNM_distance_between_samples.py -v <vcf_file> -o <output_folder> 






