# How to create a reference panel for fastNGSadmix 

## Resources
Using the following 3 files {.bed, .bim, .bed}:
http://popgen.dk/software/download/fastNGSadmix/humanOrigins_ALL/

```bash
humanOrigins_ALL/
├── humanOrigins_ALL.bed
├── humanOrigins_ALL.bim
└── humanOrigins_ALL.fam
```

## Script to generate the 2 required input files for fastNGSadmix
The  2 required input files are:
- the reference panel `refPanel`
- the file with the count of individual per population in the panel `nInd`

Steps to generate the reference panel:

```
# Download the resources  {.bed, .bim, .bed}

mkdir panel_bundle && cd panel_bundle
wget http://popgen.dk/software/download/fastNGSadmix/humanOrigins_ALL/humanOrigins_ALL.bed
wget http://popgen.dk/software/download/fastNGSadmix/humanOrigins_ALL/humanOrigins_ALL.bim
wget http://popgen.dk/software/download/fastNGSadmix/humanOrigins_ALL/humanOrigins_ALL.fam

# Inspect files:
# tree panel_bundle
#panel_bundle/
#├── humanOrigins_ALL.bed
#├── humanOrigins_ALL.bim
#└── humanOrigins_ALL.fam


# Use container with plink
docker run -it -v $PWD:$PWD -w $PWD alliecreason/plink:1.90 bash

# Install dependencies and utils
apt-get update && \
apt install r-base -y  && \
apt install  git


# Git clone fastNGSadmix repo - we will use an R script
git clone https://github.com/e-jorsboe/fastNGSadmix.git


# Convert the data into plink files:
# http://www.popgen.dk/software/index.php/FastNGSadmix

plink \
--make-bed \
--bfile humanOrigins_ALL/humanOrigins_ALL \
--out ancestry \
--maf 0.05 \
--geno 0.05 \
--chr 1-22

# output 
# $ls -l *ancestry*
# ancestry.bed
# ancestry.bim
# ancestry.fam
# ancestry.log 
```

### Installing R script dependencies
Type `R` and install from R the following packages:
- `BiocManager`
- `snpStats`

```R
install.packages("BiocManager")
BiocManager::install("snpStats")

# Exit R
quit()
```

### Convert plink files to reference panel

In the directory of the plink files, and after you have cloned the `fastNGSadmix` repo run the following command:

```bash
Rscript fastNGSadmix/R/plinkToRef.R ancestry.bed
```

You should see 3 new files generated:

```bash
ancestry.bed.sites        # 4 cols; no header; eg. 1 891021 G A
nInd_ancestry.bed.txt     # 2 rows; 1st row header (pop name); 2nd row counts of individuals per pop
refPanel_ancestry.bed.txt # ref panel we will use for both iAdmix and fastNGSadmix
```

### Post processing - Cleaning reference panel files

The files we currently use for the pipeline can be found in the following `s3` bucket and are included in `nextflow.config`



https://github.com/lifebit-ai/fast-ngs-admix/blob/828e8e2978422b471563198bc220181221addf07/nextflow.config#L1-L10

**NOTE**
The difference of these two files is the first column.
The one for iAdmix requires

```
# s3://lifebit-featured-datasets/modules/ancestry/iadmix/k27.iadmix
chr pos name A0_freq A1 AHG Amerindian Arctic ASI Australomelanesian 

# s3://lifebit-featured-datasets/modules/ancestry/fastngsadmix/refPanel_k27.txt 
id chr pos name A0_freq A1 AHG Amerindian Arctic ASI Australomelanesian
```


The generated files from the `fastNGSadmix/R/plinkToRef.R` R script, has the column `id` and many `NaN` values.

To create the reference panels with valid format, replace the `NaN` values with `0` for both `iAdmix` and `fastNGSadmix`. This heuristic might not be optimal, but for now we can test using this as a back of the envelop solution.

```bash
# replace all NaN with 0 in the reference panel
sed 's/NaN/0/g' refPanel_ancestry.bed.txt > refPanel_new.txt
```

Additionally, for generating the `iAdmix` reference panel, remove the first `id` column:

```bash
awk '{$1=""}1' refPanel_new.txt > refPanel_new_iAdmix.txt
```

The reference that has been generated now is ready.
Another optional filtering step is removing populations with very few individuals from the `nInd.txt`

```bash

#looks good
                    pop  N
Yoruba           Yoruba 70
Turkish         Turkish 56
Spanish         Spanish 53
Druze             Druze 39
Palestinian Palestinian 38
Han                 Han 33

# hmm..

                          pop N
AG2                       AG2 1
Altai                   Altai 1
Chimp                   Chimp 1
Denisova_light Denisova_light 1
Denisovan           Denisovan 1
Gorilla               Gorilla 1

```


# Newly generated reference panels in S3:

- [x] s3://lifebit-featured-datasets/modules/ancestry/iadmix/panno_no_id_col.txt
- [x] s3://lifebit-featured-datasets/modules/ancestry/fastngsadmix/panno_NaNs_replaced.txt
- [x] s3://lifebit-featured-datasets/modules/ancestry/fastngsadmix/pcount.txt


Example succesful job using those:
https://deploit.lifebit.ai/public/jobs/5dc1872adaf40500ea77688a
