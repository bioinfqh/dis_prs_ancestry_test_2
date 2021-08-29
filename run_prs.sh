#!/bin/bash

bg_dataset=$1
patient_vcf=$2
disease_list_path=$3

#patient_id=$4

ls >testoutput_new.txt

patient_id=$(python3 /scripts/get_id_from_vcf.py $patient_vcf)
bash /scripts/prs_and_report.sh $bg_dataset $patient_vcf $disease_list_path $patient_id f

