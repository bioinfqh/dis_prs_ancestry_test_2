version 1.0

task test {
  input {
    File patient_vcf
    String customer_id
    String disease_list_path
  }

  command {
  ls >testoutput_new.txt
  bash /scripts/run_prs.sh /testfiles/merged2filtered ${patient_vcf} /scripts/disease_list.txt
  cp /scripts/prs_testuser.json prs_testuser.json
  }
  output {
  File outfile = "prs_${customer_id}.json"
  }
  runtime {
  docker: "quay.io/testaccountq/dis_gen_prs_test_2:main"
  }
}

workflow make_panel_wdl {
    input {
    File patient_vcf
    }
    call test {
        input:
            patient_vcf = patient_vcf
            customer_id = "testuser",
            disease_list_path = "disease_list.txt"
            }
    output {
    test.outfile
    }
}




