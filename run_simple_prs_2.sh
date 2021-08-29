backgr_dataset=$1
dataset=$2
scores=$3
#sample_id=$4
snp_col=$4
all_col=$5
score_col=$6
prepare_dataset=$7
remove=f
filter=f
plink_loc=/plink_linux_x86_64_20201019/plink
backgr_dataset_name=${backgr_dataset##*/}

echo "scores:"
echo $scores

if [[ $scores =~ "," ]]; then
    echo "score array"
    IFS=',' read -r -a scores_array <<< "$scores"
fi
if [[ $dataset =~ ".vcf" ]]; then
    bash /scripts/convert_vcf_to_plink.sh $dataset sample1 $plink_loc sample1
    dataset_new=$(python3 /scripts/get_id_from_vcf.py $dataset)
    dataset=$dataset_new
fi

if [[ $prepare_dataset =~ "t" ]]; then
    dataset = /testfiles/sample1
    if [[ $scores =~ "," ]]; then
        touch /scripts/keep_tmp.txt
        for element in "${scores_array[@]}"
        do
        wc -l $element
        awk '{print $'$snp_col'}' $element >>/scripts/keep_tmp.txt
    done
        #wc -l keep_tmp.txt
        sort -uk1,1 /scripts/keep_tmp.txt >/scripts/keep.txt
    else
        awk '{print $'$snp_col'}' $scores >/scripts/keep.txt
    fi
    $plink_loc --bfile $dataset --extract /scripts/keep.txt --bmerge $backgr_dataset --make-bed --out /scripts/merged --allow-no-sex
    rm -f /scripts/merged.fam
    rm -f /scripts/merged.bim
    rm -f /scripts/merged.bed
    rm -f /scripts/merged-merge.fam
    rm -f /scripts/merged-merge.bim
    rm -f /scripts/merged-merge.bed
    if [ ! -f /scripts/merged-merge.missnp ]; then
        $plink_loc --bfile $dataset --extract /scripts/keep.txt --bmerge $backgr_dataset --make-bed --out /scripts/merged2 --allow-no-sex
    else
        $plink_loc --bfile $dataset --extract /scripts/keep.txt --exclude /scripts/merged-merge.missnp --make-bed --out /scripts/plink_dataset_new --allow-no-sex
        $plink_loc --bfile $backgr_dataset --extract /scripts/keep.txt --exclude /scripts/merged-merge.missnp --make-bed --out /scripts/data_file_new --allow-no-sex
        $plink_loc --bfile /scripts/plink_dataset_new --extract /scripts/keep.txt --bmerge /scripts/data_file_new --make-bed --out /scripts/merged --allow-no-sex
        rm -f /scripts/data_file_new.bed
        rm -f /scripts/data_file_new.bim
        rm -f /scripts/data_file_new.fam
        rm -f /scripts/plink_dataset_new.bed
        rm -f /scripts/plink_dataset_new.bim
        rm -f /scripts/plink_dataset_new.fam
        $plink_loc --bfile /scripts/merged --exclude /scripts/merged-merge.missnp --make-bed --out /scripts/merged2 --allow-no-sex
    fi
    rm -f /scripts/merged.fam
    rm -f /scripts/merged.bim
    rm -f /scripts/merged.bed
    rm -f /scripts/merged-merge.fam
    rm -f /scripts/merged-merge.bim
    rm -f /scripts/merged-merge.bed
    bfile_loc=/scripts/merged2
else
    bfile_loc=$backgr_dataset
    bfile_loc=/testfiles/testsample
fi

if [[ $scores =~ "," ]]; then
    ctr=1
    for element in "${scores_array[@]}"
    do
        echo "$element"
        $plink_loc --bfile $bfile_loc --score $element $snp_col $all_col $score_col --out /scripts/scores_$ctr
        #$plink_loc --bfile merged2 --score scores_new_$ctr $snp_col $all_col $score_col --extract outfile_tmp.valid.snp --out outfile_tmp_$ctr
        #sudo rm -f summary_stats_filtered$ctr.pvalue
        let ctr++
    done
else
    echo $scores
    #$plink_loc --bfile merged2 --score $scores $snp_col $all_col $score_col --extract outfile_tmp.valid.snp --out outfile_tmp
    $plink_loc --bfile $bfile_loc --score $scores $snp_col $all_col $score_col --out /scripts/scores_1
fi
if [[ $remove == "t" ]]; then
    rm -f /scripts/merged2.fam
    rm -f /scripts/merged2.bim
    rm -f /scripts/merged2.bed
    rm -f /scripts/outfile_tmp.valid.snp
fi
rm -f /scripts/merged-merge.missnp
rm -f /scripts/summary_stats_filtered.pvalue
rm -f /scripts/keep.txt
rm -f /scripts/clumped.clumped

result_id=${dataset##*/}
#if [[ $scores =~ "," ]]; then
#    ctr=0
#    for element in "${scores_array[@]}"
#    do
#        python3 read_results.py scores_$ctr.0.1.profile $result_id 2 6
#        let ctr++
#    done
#else
#    python3 read_results.py scores_test.0.1.profile $result_id 2 6
#fi


#$plink_loc --bfile $dataset --score $scores $snp_col $all_col $score_col --extract outfile_tmp.valid.snp --out outfile_tmp

#python3 read_results.py /home/quirin/Downloads/outfile_test_99.profile 3285b0mtsjwbtxgwlukasqhflke3_fc3lh9p 2 6
