
dataset=$1
summary_stats=$2
snp_col=$3
all_col=$4
score_col=$5
pval_col=$6
use_log=$7
plink_loc=/home/quirin/Downloads/plink_linux_x86_64_20201019/plink 
backgr_dataset=/media/quirin/INTENSO/merged1K
#tmpdir_path=tmp_dir

if [[ $dataset =~ ".vcf" ]]; then
    bash convert_vcf_to_plink.sh $dataset input_file_tmp $plink_loc sample_1
    dataset=input_file_tmp
fi

if [[ $use_log =~ "log" ]]; then
    python3 logarithmize.py $summary_stats $score_col >summary_stat_tmp.txt
    summary_stats=summary_stat_tmp.txt
fi

head -n 1 $summary_stats >firstline.txt
awk '{$'$pval_col' = "P"; print}' firstline.txt >summary_stats_tmp_2.txt
awk '{$'$snp_col' = "SNP"; print}' summary_stats_tmp_2.txt >summary_stats_converted.txt
tail -n +1 $summary_stats >>summary_stats_converted.txt

$plink_loc --bfile $backgr_dataset --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump summary_stats_converted.txt --clump-snp-field SNP --clump-field P --out /media/quirin/INTENSO/clumped

awk 'NR!=1{print $'$snp_col'}' /media/quirin/INTENSO/clumped.clumped >outfile_tmp.valid.snp
#awk '!seen[$'$snp_col']++' summary_stats_converted.txt >summary_stats_filtered.txt
sort -uk$snp_col,$snp_col summary_stats_converted.txt >summary_stats_filtered.txt

awk '{print $'$snp_col'}' summary_stats_filtered.txt >keep.txt

$plink_loc --bfile $backgr_dataset --extract keep.txt --make-bed --out filtered_$dataset --allow-no-sex
awk '{print $'$snp_col',$'$pval_col'}' summary_stats_filtered.txt > summary_stats_filtered.pvalue

echo "0.001 0 0.001" > range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

#$plink_loc --bfile filtered_$dataset --score summary_stats_filtered.txt $snp_col $all_col $score_col --q-score-range range_list summary_stats_filtered.pvalue --extract outfile_tmp.valid.snp --out scores_$dataset

$plink_loc --bfile $dataset --extract keep.txt --bmerge filtered_$dataset --make-bed --out merged --allow-no-sex
sudo rm -f merged.fam
sudo rm -f merged.bim
sudo rm -f merged.bed
sudo rm -f merged-merge.fam
sudo rm -f merged-merge.bim
sudo rm -f merged-merge.bed
$plink_loc --bfile $dataset --extract keep.txt --exclude merged-merge.missnp --make-bed --out plink_dataset_new --allow-no-sex
$plink_loc --bfile filtered_$dataset --extract keep.txt --exclude merged-merge.missnp --make-bed --out data_file_new --allow-no-sex
$plink_loc --bfile plink_dataset_new --extract keep.txt --bmerge data_file_new --make-bed --out merged --allow-no-sex
rm -f data_file_new.bed
rm -f data_file_new.bim
rm -f data_file_new.fam
rm -f plink_dataset_new.bed
rm -f plink_dataset_new.bim
rm -f plink_dataset_new.fam
sudo $plink_loc --bfile merged --exclude merged-merge.missnp --make-bed --out merged2 --allow-no-sex
sudo rm -f merged.fam
sudo rm -f merged.bim
sudo rm -f merged.bed
sudo rm -f merged-merge.fam
sudo rm -f merged-merge.bim
sudo rm -f merged-merge.bed
$plink_loc --bfile merged2 --score summary_stats_filtered.txt $snp_col $all_col $score_col --q-score-range range_list summary_stats_filtered.pvalue --extract outfile_tmp.valid.snp --out scores_final_$dataset
#$plink_loc --bfile merged2 --score summary_stats_filtered.txt $snp_col $all_col $score_col --extract outfile_tmp.valid.snp --out outfile_tmp
sudo rm -f merged2.fam
sudo rm -f merged2.bim
sudo rm -f merged2.bed
sudo rm -f filtered_$dataset.bed
sudo rm -f filtered_$dataset.bim
sudo rm -f filtered_$dataset.fam
sudo rm -f outfile_tmp.valid.snp
sudo rm -f merged-merge.missnp
sudo rm -f summary_stats_filtered.pvalue
sudo rm -f keep.txt
sudo rm -f summary_stats_filtered.txt
sudo rm -f /media/quirin/INTENSO/clumped.clumped

result_id=${dataset##*/}
python3 read_results.py scores_test.0.1.profile $result_id 2 6
