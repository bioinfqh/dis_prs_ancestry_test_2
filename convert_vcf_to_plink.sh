
input_path=$1
output_path=$2
plink_path="/plink_linux_x86_64_20201019/plink"
name=${input_path##*/}
name=$4


name=$(python3 /scripts/get_id_from_vcf.py $input_path)
$plink_path --vcf $input_path --make-bed --out /scripts/$name
echo NA $name 0 0 0 1 >/scripts/$name.fam
