version 1.0

task test {
  input {
    String patient_vcf
    String customer_id
    String disease_list_path
  }

  command {
  ls >output_new.txt
  }
  output {
  File outfile = "output_new.txt"
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

