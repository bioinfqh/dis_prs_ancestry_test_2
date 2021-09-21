
patient_vcf=$1
disease_list_path=$2
customer_id=$3


outfiles=$(python3 /scripts/extract_table.py $patient_vcf all /scripts/dis_report_$customer_id $customer_id)

echo $outfiles

touch /scripts/list_of_outfiles_$customer_id.txt

json_posted="false"

json_posted=$(python3 /scripts/post_json.py $outfiles multiple_files)
echo $json_posted
#IFS='_SEPARATOR_' read -ra ADDR <<< "$outfiles"
#for i in "${ADDR[@]}"; do
#  echo $i
#  json_posted=$(python3 /scripts/post_json.py $i)
#  echo $i >>/scripts/list_of_outfiles_$customer_id.txt
#  echo $json_posted
#done

