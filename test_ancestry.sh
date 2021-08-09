#!/bin/bash



## check if settings from settings.txt are used
if [[ -f settings.txt ]]; then
anc_line=$(grep ancestry settings.txt)
ancestry_path=${anc_line##*=}
refpanel_line=$(grep refpanel settings.txt)
refpanel_path=${refpanel_line##*=}
plink_line=$(grep plink settings.txt)
plink_path=${plink_line##*=}
elif [[ $plink_path == "" ]];then
ancestry_path="ancestry_mirror-master/src/ancestry"
refpanel_path="reference_files/refpanel_26pops.txt.gz"
plink_path="plink"
fi


## check for show plot option
genotype_file=dis_calc_application_cont/testfiles/file_for_prs.vcf
outfile=outfile_testuser
outfile_results=outfile_results
params=""

# define temporary file paths
out="output"
tmp=${genotype_file##*/}
dtc_file=${tmp//.txt/}
gzip_tmp="gzip_tmp.txt.gz"
this_dir=$PWD
#refpanel="reference_files/refpanel_26pops.txt.gz"


## convert files to gzipped 23andMe file
if [[ $genotype_file =~ zip ]]; then
      unzip $genotype_file
      gzip -cdf $genotype_file > $genotype_file
fi
if [[ $genotype_file =~ "vcf" ]]; then
    $plink_path --vcf $genotype_file --double-id --snps-only --biallelic-only --recode 23 --out out_dir 
    #> /dev/null 2>&1
    gzip <out_dir.txt >out_dir/gzip_tmp.txt.gz
elif [[ $genotype_file =~ ".23andMe.txt" ]]; then
    gzip -cdfq $genotype_file > out_dir/$gzip_tmp
elif [[ $genotype_file =~ "txt" ]]; then
    python3 docker/common-latest-geno_ancestry/dtc_genomics_convertor.py -i $genotype_file -o out_dir/$dtc_file -f txt  > /dev/null 2>&1
    gzip -cvfq out_dir/$dtc_file.txt  > out_dir/$gzip_tmp
fi

export LD_LIBRARY_PATH=/htslib-1.3.2

#echo $params
## run ancestr
#params2="$params"

$ancestry_path -i $refpanel_path -g23 out_dir/$gzip_tmp -e 0.01 -o output $params 
if [[ -f output.Q ]]; then
    echo "output written to output.Q"
else
    echo "output not written to file"
fi
#> /dev/null 2>&1
#$ancestry_path -i $refpanel_path -g23 out_dir/$gzip_tmp "${params}" -o $out 
#$ancestry_path -i $refpanel_path -g23 out_dir/$gzip_tmp -e 0.01 --merge -o $out > /dev/null 2>&1


## delete temporary files
if [ -f out_dir/$gzip_tmp ]; then
rm out_dir/$gzip_tmp
fi

if [ -f out_dir/$dtc_file ]; then
rm out_dir/$dtc_file
fi
rm -r -f work

id=${dtc_file//.23andMe.txt/}

### this prints: ['italian': 0.5, 'spanish': 0.5 ] }
printf "[ "
while read p; do
  printf "'"
  echo -n $p | awk '{print $1}'
  printf "':"
  echo -n $p | awk '{print $2}' 
  printf " , " 
done<output.Q 
printf " ]"

Rscript -e 'ancestry <- read.table( "output.Q" , sep = " ");  data.table::fwrite(as.data.frame(t(ancestry)), sep = ",", quote = FALSE, file = "ancestry.csv", col.names = FALSE)'

if [[ -f ancestry.csv ]]; then
    echo "ancestry table written to ancestry.csv"
else
    echo "ancestry table not written to ancestry.csv"
fi

# copy the docker bin into pwd
#mkdir bin/
#cp -r /opt/conda/envs/fast-ngs-admix/bin/* bin/
cd bin

  # copy the rmarkdown into the workdir
R -e "rmarkdown::render('map.Rmd', params = list(ancestry_csv='../ancestry.csv'), output_file='map.html')"

if [[ -f map.html ]]; then
    echo "map.html was created."
else
    echo "map.html was not created."
fi


cd ..



mkdir MultiQC
cp bin/map.html MultiQC/multiqc_report.html
cp output.Q $outfile_results
mv bin/map.html $outfile
#mv bin/map.html dis_calc/static/ancestry_map.html
### this prints only the ancestry proportions, line wise
#while read p; do
#  echo $p
#done<output.Q

