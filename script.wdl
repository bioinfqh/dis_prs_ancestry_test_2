version 1.0

task test {
  input {
    File patient_vcf
    String customer_id
    String disease_list_path
  }

  command {
  bash dis_calc/run_prs.sh patient_vcf patient_vcf disease_list_path
  }
  output {
  File outfile = "dis_calc/static/prs_${customer_id}.json"
  }
  runtime {
  docker: "https://github.com/bioinfqh/dis_prs_ancestry"
  }
}

workflow make_panel_wdl {
    call test {
        input:
            patient_vcf = prs_vcf.vcf,
            customer_id = "testuser",
            disease_list_path = "disease_list.txt"
            }
}

