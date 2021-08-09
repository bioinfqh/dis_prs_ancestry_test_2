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
plink_loc=plink_linux_x86_64_20201019/plink
backgr_dataset_name=${backgr_dataset##*/}

echo "scores:"
echo $scores

if [[ $scores =~ "," ]]; then
    echo "score array"
    IFS=',' read -r -a scores_array <<< "$scores"
fi
if [[ $dataset =~ ".vcf" ]]; then
    bash convert_vcf_to_plink.sh $dataset sample1 $plink_loc sample1
    dataset_new=$(python3 get_id_from_vcf.py $dataset)
    dataset=$dataset_new
fi

if [[ $prepare_dataset =~ "t" ]]; then

    if [[ $scores =~ "," ]]; then
        touch keep_tmp.txt
        for element in "${scores_array[@]}"
        do
        wc -l $element
        awk '{print $'$snp_col'}' $element >>keep_tmp.txt
    done
        #wc -l keep_tmp.txt
        sort -uk1,1 keep_tmp.txt >keep.txt
    else
        awk '{print $'$snp_col'}' $scores >keep.txt
    fi
    $plink_loc --bfile $dataset --extract keep.txt --bmerge $backgr_dataset --make-bed --out merged --allow-no-sex
    rm -f merged.fam
    rm -f merged.bim
    rm -f merged.bed
    rm -f merged-merge.fam
    rm -f merged-merge.bim
    rm -f merged-merge.bed
    if [ ! -f merged-merge.missnp ]; then
        $plink_loc --bfile $dataset --extract keep.txt --bmerge $backgr_dataset --make-bed --out merged2 --allow-no-sex
    else
        $plink_loc --bfile $dataset --extract keep.txt --exclude merged-merge.missnp --make-bed --out plink_dataset_new --allow-no-sex
        $plink_loc --bfile $backgr_dataset --extract keep.txt --exclude merged-merge.missnp --make-bed --out data_file_new --allow-no-sex
        $plink_loc --bfile plink_dataset_new --extract keep.txt --bmerge data_file_new --make-bed --out merged --allow-no-sex
        rm -f data_file_new.bed
        rm -f data_file_new.bim
        rm -f data_file_new.fam
        rm -f plink_dataset_new.bed
        rm -f plink_dataset_new.bim
        rm -f plink_dataset_new.fam
        $plink_loc --bfile merged --exclude merged-merge.missnp --make-bed --out merged2 --allow-no-sex
    fi
    rm -f merged.fam
    rm -f merged.bim
    rm -f merged.bed
    rm -f merged-merge.fam
    rm -f merged-merge.bim
    rm -f merged-merge.bed
    bfile_loc=merged2
else
    bfile_loc=$backgr_dataset
fi

if [[ $scores =~ "," ]]; then
    ctr=1
    for element in "${scores_array[@]}"
    do
        echo "$element"
        $plink_loc --bfile $bfile_loc --score $element $snp_col $all_col $score_col --out scores_$ctr
        #$plink_loc --bfile merged2 --score scores_new_$ctr $snp_col $all_col $score_col --extract outfile_tmp.valid.snp --out outfile_tmp_$ctr
        #sudo rm -f summary_stats_filtered$ctr.pvalue
        let ctr++
    done
else
    echo $scores
    #$plink_loc --bfile merged2 --score $scores $snp_col $all_col $score_col --extract outfile_tmp.valid.snp --out outfile_tmp
    $plink_loc --bfile $bfile_loc --score $scores $snp_col $all_col $score_col --out scores_1
fi
if [[ $remove == "t" ]]; then
    rm -f merged2.fam
    rm -f merged2.bim
    rm -f merged2.bed
    rm -f outfile_tmp.valid.snp
fi
rm -f merged-merge.missnp
rm -f summary_stats_filtered.pvalue
rm -f keep.txt
rm -f /media/quirin/INTENSO/clumped.clumped

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
