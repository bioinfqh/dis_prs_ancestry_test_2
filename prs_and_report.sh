
backgr_dataset=$1
dataset=$2
#param_file=$3
disease_list_file=$3
patient_id=$4
predict=$5


# predict is set to true if absolute risk should be calculated 
# if backgr_dataset and dataset are identical, the script assumes that the patient is already in the reference file and uses this file.
ctr=1
rm -f /scripts/resultparamfile.txt
touch /scripts/resultparamfile.txt
phenos=()
#scores=()
disease_names=()
### read disease names, paths of SNP wise score files, path of phenotype files
while read line; do
    disease_name=$(echo -n $line | awk -F',' '{print $1}')
    #echo $disease_name
    score_curr=$(echo -n $line | awk -F',' '{print $2}')
    pheno_curr=$(echo -n $line | awk -F',' '{print $3}')
    #tmp=$disease_name,
    if [[ $ctr -eq 1 ]]; then
        disease_names+=($disease_name)
        scores=$score_curr
        #scores+=($score_curr)
        phenos+=($pheno_curr)
    else
        disease_names+=($disease_name)
        #scores+=($score_curr)
        scores="$scores,$score_curr"
        phenos+=($pheno_curr)
        #phenos="${phenos},$pheno_curr"
    fi
    # write file with paths to PRS result files
    echo -e $disease_name'\t'/scripts/scores_$ctr.profile >>/scripts/resultparamfile.txt
    let ctr++
done<$disease_list_file

echo "giving scores"
echo $scores
#sudo bash run_simpe_prs_2.sh $backgr_dataset $dataset $scores 1 2 4 t
## run PRS calculation
if [[ $dataset == $backgr_dataset ]]; then
    bash /scripts/run_simple_prs_2.sh $backgr_dataset $dataset $scores 1 2 3 f
else
    bash /scripts/run_simple_prs_2.sh $backgr_dataset $dataset $scores 1 2 3 t
fi
# get patient ID from VCF
patient_id_new=$(python3 /scripts/get_id_from_vcf.py $dataset)

#bash dis_calc/run_simple_prs_2.sh $backgr_dataset $dataset $scores 1 2 3 3 
rm -f /scripts/results_for_script.txt
touch /scripts/results_for_script.txt
ctr=1
#echo "${phenos[*]}"
#echo ${phenos[0]}
#echo ${phenos[1]}
for dis in "${disease_names[@]}"
do
    if [[ $predict == "f" ]]; then
        # calculate relative risk and write to result file
        tmp=$(python3 /scripts/get_relative_risk.py /scripts/scores_$ctr.profile $patient_id_new /scripts/results_and_phenos$ctr.txt)
        risk=$(echo $tmp | awk '{print $2}')
        percentile=$(echo $tmp | awk '{print $3}')
    else
        # calculate absolute risk based on phenotype data of reference samples
        tmp=$(bash /scripts/get_absolute_risk.sh /scripts/scores_$ctr.profile ${phenos[ctr]} $patient_id_new /scripts/results_and_phenos$ctr.txt)
        risk=$(echo $tmp | awk '{print $2}')
        percentile=$(echo $tmp | awk '{print $3}')
    fi
    echo -e $dis"\t"$risk"\t"$percentile"\t"/scripts/results_and_phenos$ctr.txt >>/scripts/results_for_script.txt
    let ctr++
done

## write report
python3 /scripts/make_prs_report.py "easy" /scripts/results_for_script.txt $patient_id_new "your_prs_report.pdf" $predict
#sudo python3 make_prs_report.py resultparamfile.txt $patient_id "your_prs_report.pdf"
#ctr2=0
#for (( c=0; c<=$ctr; c++ ))
#do
#    echo 
#done


