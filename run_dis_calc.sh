
patient_vcf=$1
disease_list_path=$2
customer_id=$3

outfiles=$(python3 /scripts/extract_table.py $patient_vcf $disease_list_path /scripts/dis_report_$customer_id)

touch /scripts/list_of_outfiles_$customer_id.txt

IFS='_SEPARATOR_' read -ra ADDR <<< "$outfiles"
for i in "${ADDR[@]}"; do
  echo $i
  echo $i >>/scripts/list_of_outfiles_$customer_id.txt
done
