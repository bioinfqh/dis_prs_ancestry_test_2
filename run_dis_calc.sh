
patient_vcf=$1
disease_list_path=$2
customer_id=$3

outfiles=$(python3 extract_table.py $patient_vcf $disease_list_path dis_calc/static/dis_report_$customer_id)

touch list_of_outfiles_$customer_id.txt

IFS='_SEPARATOR_' read -ra ADDR <<< "$outfiles"
for i in "${ADDR[@]}"; do
  echo $i
  echo $i >>list_of_outfiles_$customer_id.txt
done
