version 1.0

task test {
  input {
    String patient_vcf
    String customer_id
    String disease_list_path
  }

  command {
  ls >testoutput_new.txt
  bash /scripts/run_prs.sh patient_vcf patient_vcf disease_list_path
  }
  output {
  File outfile = "prs_${customer_id}.json"
  }
  runtime {
  docker: "quay.io/testaccountq/dis_gen_prs_test_2:main"
  }
}

workflow make_panel_wdl {
    call test {
        input:
            patient_vcf = "prs_vcf.vcf",
            customer_id = "testuser",
            disease_list_path = "disease_list.txt"
            }
}



