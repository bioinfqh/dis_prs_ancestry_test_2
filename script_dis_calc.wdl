version 1.0

task test {
  input {
    File patient_vcf
    File? file_with_varsome_terms = "/scripts/list_of_terms_for_varsome.txt"
    String disease_list
    String customer_id
  }

  command {
  ls >testoutput_new.txt
  cp ${file_with_varsome_terms} /scripts/list_of_terms_for_varsome_2.txt
  mv /scripts/list_of_terms_for_varsome_2.txt /scripts/list_of_terms_for_varsome.txt
  bash /scripts/run_dis_calc.sh ${patient_vcf} /scripts/disease_groups_dis_calc.txt testuser
  cp /scripts/dis_genes_testuser_all.json dis_genes_testuser_all.json
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
    File? file_with_varsome_terms = "/scripts/list_of_terms_for_varsome.txt"
    }
    call test {
        input:
            patient_vcf = patient_vcf,
            customer_id = "testuser",
            disease_list = "disease_list.txt"
            file_with_varsome_terms = file_with_varsome_terms
            }
    output {
    File resultfile = test.outfile
    }
}


