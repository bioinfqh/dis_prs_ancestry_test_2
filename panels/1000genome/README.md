# 26 populations from 1000genome

The function [generatePopLists.R](https://github.com/cgpu/fast-ngs-admix/blob/master/bin/generatePopLists.R) can be used to create lists of the sample names that belond in each subpopulation.


# Metadata table from 1000genomes

The table used for generating the lists can be found in the following link:

```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```
And is also uploaded with public read access here on our s3:

```
s3://fast-ngs/integrated_call_samples_v3.20130502.ALL.panel
```

The table looks like this:

```
sample  pop     super_pop       gender
HG00096 GBR     EUR     male
HG00097 GBR     EUR     female
.
.
```


# S3 Bucket (public read access)

For the 2504 from 1000genome, the generated lists with sample names are store in this folder:
```
s3://fast-ngs/26subpopulations/
```

```
├── african_american_sw.txt
├── african_caribbean.txt
├── bengali.txt
├── british.txt
├── ceph.txt
├── colombian.txt
├── dai_chinese.txt
├── esan.txt
├── finnish.txt
├── gambian.txt
├── gujarati.txt
├── han_chinese.txt
├── indian.txt
├── japanese.txt
├── kinh_vietnamese.txt
├── luhya.txt
├── mende.txt
├── mexican_american.txt
├── peruvian.txt
├── puerto_rican.txt
├── punjabi.txt
├── southern_han_chinese.txt
├── spanish.txt
├── sri_lankan.txt
├── tuscan.txt
└── yoruba.txt
```