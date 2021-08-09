# `lists`

> Sample names (eg HG00096) for each one of the 26 populations of 1000genomes

The folder [lists](https://github.com/lifebit-ai/fast-ngs-admix/tree/master/panels/1000genome/populations/lists) contains single column .txt files with the sample ids that correspond to each one of the 26 populations from the 1000genomes project. 

The files look like this:

```bash
cd lists && head -2 *txt
# ==> african_american_sw_n_61.txt <==
# NA19625
# NA20359

# ==> african_caribbean_n_96.txt <==
# HG01914
# HG02144
```


_Q: How can I generate those lists?_

A: With the function [generatePopLists.R](https://github.com/cgpu/fast-ngs-admix/blob/master/bin/generatePopLists.R)

```R
library(data.table)
source("generatePopLists.R")

pop_metadata <- data.table::fread("metadata_populations_1000genomes.csv")

# > head(pop_metadata, 2)
#    pop  sample super_pop gender   population_name
# 1: ACB HG01914       AFR   male african_caribbean
# 2: ACB HG02144       AFR female african_caribbean

# This will populate `savedir` with .txt files with sample ids
generatePopLists( whole_pop_df = augmented_country_pop_codes_df, 
                  sample_col   = "sample", 
                  pop_col      = "population_name",  
                  savedir      = "~/some/where")
```

The populated folder is similar to [lists](https://github.com/lifebit-ai/fast-ngs-admix/tree/master/panels/1000genome/populations/lists), depending on the input metadata `.csv`.
The counts indicate the number of individuals in each population.

```bash
cd ~/some/where
tree 
# ├── african_american_sw_n_61.txt
# ├── african_caribbean_n_96.txt
# ├── bengali_n_86.txt
# ├── british_n_91.txt
# ├── ceph_n_99.txt
# ├── colombian_n_94.txt
# ├── dai_chinese_n_93.txt
# ├── esan_n_99.txt
# ├── finnish_n_99.txt
# ├── gambian_n_113.txt
# ├── gujarati_n_103.txt
# ├── han_chinese_n_103.txt
# ├── indian_n_102.txt
# ├── japanese_n_104.txt
# ├── kinh_vietnamese_n_99.txt
# ├── luhya_n_99.txt
# ├── mende_n_85.txt
# ├── mexican_american_n_64.txt
# ├── peruvian_n_85.txt
# ├── puerto_rican_n_104.txt
# ├── punjabi_n_96.txt
# ├── southern_han_chinese_n_105.txt
# ├── spanish_n_107.txt
# ├── sri_lankan_n_102.txt
# ├── tuscan_n_107.txt
# └── yoruba_n_108.txt
```


# `vcfs`

> Subsetting the 1000genomes multisample VCF files

The sample IDs stored in the `lists` folder, are needed to subset the VCF files from 1000genomes, which are provided per chromosome and include the whole cohort of N = 2504 individuals.
The VCF files are available here:

```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
```

In the folder `vcfs` there is a `.csv` file tha includes all autosomal chromosomes (2 columns, cols ["vcf","index"]; bgzipped files with `vcf.gz` suffix, and respective indices `vcf.gz.tbi`).  We choose to exclude sex chromosome snps from the ancestry reference panel.

The table is a required input for [`lifebit-ai/panelR`](https://github.com/lifebit-ai/panelR) (`--multiVCF_table  vcfs/1000genomes_autosomal_chr_vcfs.csv` ).
Here is a preview of how the file looks like:

```bash
head -2  vcfs/1000genomes_autosomal_chr_vcfs.csv

# added {..} to shorten path for docs preview, real path in the csv file
# vcf,index
# ftp://ftp.1000genomes.ebi.ac.uk/{..}.chr1.{..}.vcf.gz,ftp://ftp.1000genomes.ebi.ac.uk/{..}.chr1.{..}.vcf.gz.tbi
# ftp://ftp.1000genomes.ebi.ac.uk/{..}.chr2.{..}.vcf.gz,ftp://ftp.1000genomes.ebi.ac.uk/{..}.chr2.{..}.vcf.gz.tbi
```


