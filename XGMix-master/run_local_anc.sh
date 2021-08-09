infile=$1
sample_id=$2
customer_id=$3

infile_dir="/media/quirin/INTENSO/XGMIX_infiles"
model_file_path="/media/quirin/INTENSO/XGMIX_model_files"
output_prefix_path="demo_data/demo"

#bash partition_vcf.sh $infile vcf_temp_dir
outfile=$(python3 run_xgmix.py vcf_temp_dir $sample_id $model_file_path $output_prefix_path $customer_id)
