version 1.0

task test {
  input {
    File patient_vcf
    String disease_list
    String customer_id
  }

  command {
  ls >testoutput_new.txt
  bash /scripts/run_dis_calc.sh ${patient_vcf} /scripts/disease_groups_dis_calc.txt ${customer_id}
  cp /scripts/dis_genes_${customer_id}_all.json dis_genes_${customer_id}_all.json
  }
  output {
  File outfile = "dis_genes_${customer_id}_all.json"
  }
  runtime {
  docker: "quay.io/testaccountq/dis_gen_prs_test_2:main"
  }
}

workflow make_panel_wdl {
    input {
    File patient_vcf
    String? customer_id = "testuser"
    }
    call test {
        input:
            patient_vcf = patient_vcf,
            customer_id = customer_id,
            disease_list = "disease_list.txt"
            }
    output {
    File resultfile = test.outfile
    }
}

