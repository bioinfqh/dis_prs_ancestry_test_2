version 1.0

task test {
  input {
    String input_vcf
    String customer_id
  }

  command {
  ls >testoutput_new.txt
  cp input_vcf /scripts/anc_vcf.vcf
  bash /scripts/run_global_anc.sh /testfiles/file_for_prs.vcf ${customer_id}
  bash /scripts/run_local_anc.sh /testfile/file_for_prs.vcf test ${customer_id}
  cp /scripts/ancestry_testuser.json ancestry_testuser.json
  cp /scripts/chm_img_testuser.png chm_img_testuser.png
  }
  output {
  File outfile = "/scripts/ancestry_${customer_id}.json"
  File outfile2 = "/scripts/chm_img_${customer_id}.png"
  }
  runtime {
  docker: "quay.io/testaccountq/dis_gen_prs_test_2:main"
  }
}

workflow make_panel_wdl {
    call test {
        input:
            input_vcf = "prs_vcf.vcf",
            customer_id = "testuser"
            }
}




