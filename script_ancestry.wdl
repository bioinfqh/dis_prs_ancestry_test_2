version 1.0

task test {
  input {
    File patient_vcf
    String customer_id
  }

  command {
  ls >testoutput_new.txt
  cp ${patient_vcf} /scripts/anc_vcf.vcf
  bash /scripts/run_global_anc.sh /scripts/anc_vcf.vcf ${customer_id}
  bash /scripts/run_local_anc.sh /scripts/anc_vcf.vcf test ${customer_id}
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
input {
    File patient_vcf
    }
    call test {
        input:
            patient_vcf = patient_vcf,
            customer_id = "testuser"
            }
    output {
    File resultfile = test.outfile
    }
}
