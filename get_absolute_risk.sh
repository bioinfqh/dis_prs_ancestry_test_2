
gwas_results=$1
phenos=$2
patient_id=$3
outfile=$4

sort -k2 $gwas_results >gwas_results_sorted.txt
sort -k1 $phenos >phenos_sorted.txt
join -1 2 -2 1 -o 1.1 1.6 2.2 gwas_results_sorted.txt phenos_sorted.txt >prs_and_pheno.txt
sort -k2 prs_and_pheno.txt >results_temp.txt
pred_line=$(python3 risk_calc_2.py results_temp.txt $patient_id $outfile)
#mv results_temp.txt $outfile
echo $pred_line
#abs_risk = $(echo -n $pred_line | awk -F',' '{print $2}')
#percentile = $(echo -n $pred_line | awk -F',' '{print $3}')
