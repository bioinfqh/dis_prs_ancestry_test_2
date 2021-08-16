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
  docker: "github.com/bioinfqh/dis_prs_ancestry_test_2"
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

