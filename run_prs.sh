#!/bin/bash

bg_dataset=$1
patient_vcf=$2
disease_list_path=$3

#patient_id=$4

patient_id=$(python3 get_id_from_vcf.py $patient_vcf)
bash prs_and_report.sh $bg_dataset $disease_list_path $patient_id f

