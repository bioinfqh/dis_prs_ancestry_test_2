infile=$1
sample_id=$2
customer_id=$3

infile_dir="/XGMIX_infiles"
model_file_path="/XGMIX_model_files"
output_prefix_path="/XGMix_master/demo_data/demo"

bash /scripts/partition_vcf.sh $infile /vcf_temp_dir
outfile=$(python3 /scripts/run_xgmix.py /vcf_temp_dir $sample_id $model_file_path $output_prefix_path $customer_id)
