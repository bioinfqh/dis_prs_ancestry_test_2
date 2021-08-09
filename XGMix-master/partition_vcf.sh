input_vcf=$1
outdir=$2

bgzip -c $input_vcf > myvcf.vcf.gz
tabix -p vcf myvcf.vcf.gz

for i in 21 22
do
    bcftools view -r $i myvcf.vcf.gz >$outdir/chr_$i.vcf
    #tabix myvcf.vcf.gz chr$i > test_dir/chr$i.vcf
done

rm myvcf.vcf.gz
