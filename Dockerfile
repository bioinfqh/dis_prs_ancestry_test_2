FROM ubuntu:18.04
COPY 23dir 23dir
COPY ancestry_mirror-master ancestry_mirror-master
COPY assets assets
COPY converted converted
COPY example example
COPY gsl-2.6 gsl-2.6
#COPY MultiQC MultiQC
COPY out_dir out_dir
COPY output output
COPY output.23andMe output.23andMe
COPY panels panels
COPY pops pops
COPY reference_files reference_files
COPY testdata testdata
COPY testfiles testfiles
COPY bin bin
#COPY XGMIX_infiles XGMIX_infiles
COPY XGMix-master XGMix-master
COPY XGMIX_model_files XGMIX_model_files
COPY plink_linux_x86_64_20201019 plink_linux_x86_64_20201019
COPY resultparamfile.txt .
COPY results_for_script.txt .
COPY syndict_temp.txt .
COPY urls.py .
COPY get_ancestry.sh .
COPY convert.py .
COPY convert_nats.py .
COPY install.sh .
COPY install_2.sh .
COPY install_packages.R .
COPY setup.py .
COPY summarize.py .
COPY settings.txt .
COPY gsl-2.6.tar.gz .
COPY htslib.tar.bz2 .
COPY install_packages_2.sh .
COPY install_packages_2.R .
COPY install_tidyverse.R .
COPY requirements.txt .
COPY base.svg.p .
COPY run_local_anc.sh .
COPY script.wdl .
COPY test.json .
COPY extract_table.py .
COPY get_absolute_risk.sh .
COPY get_id_from_vcf.py .
COPY get_relative_risk.py .
COPY list_of_prs.txt .
COPY logarithmize.py .
COPY make_prs_report.py .
COPY make_report_from_files.py .
COPY partition_vcf.sh .
COPY plot_curve.py .
COPY pop_list.txt .
COPY prepare_plots.py .
COPY prs_and_report.sh .
COPY read_results.py .
COPY rsid_to_clinvar.py .
COPY run_dis_calc.sh .
COPY run_docker.sh .
COPY run_local_anc.sh .
COPY run_prs.sh .
COPY run_simple_prs_2.sh .
COPY run_xgmix.py .
COPY setup.py . 
COPY simple_prs.sh .
COPY tableExport.csv .
COPY test_ancestry.sh .
COPY tests.py .
COPY urls.py .
COPY views.py .
COPY views_2.py .
#RUN apt-get update
RUN mkdir MultiQC
RUN mkdir XGMIX_infiles
RUN chmod 777 plink_linux_x86_64_20201019/* && chmod 777 ancestry_mirror-master/src/ancestry 
#RUN apt-get update && apt-get install -y  python3.7 && apt-get install -y python3-pip && apt autoremove && apt autoclean && apt-get clean packages && pip3 install -r requirements.txt && apt autoremove && apt autoclean && apt-get clean packages && apt-get install -y curl && python3 manage.py makemigrations && python3 manage.py migrate && python3 manage.py runserver 0.0.0.0:8000
#RUN apt-get update && apt-get install -y  python3.7 && apt-get install -y python3-pip && apt autoremove && apt autoclean && apt-get clean packages && pip3 install --upgrade setuptools && pip3 install -r requirements.txt && apt autoremove && apt autoclean && apt-get clean packages && apt-get install -y curl && apt-get install -y wget && wget https://github.com/wkhtmltopdf/wkhtmltopdf/releases/download/0.12.4/wkhtmltox-0.12.4_linux-generic-amd64.tar.xz && tar -xvf wkhtmltox-0.12.4_linux-generic-amd64.tar.xz && mv wkhtmltox/bin/wkhtmltopdf /usr/bin/wkhtmltopdf &&  mv wkhtmltox/bin/wkhtmltoimage /usr/bin/wkhtmltoimage && apt-get install -y libxrender1 && apt-get install -y libfontconfig && apt-get install -y libxext6 && apt-get install -y libssl1.0-dev && DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && apt-get install -y r-base && apt-get install -y sudo && apt-get install -y make && bash install.sh && tar -xjvf htslib.tar.bz2 && export LD_LIBRARY_PATH=/usr/local/lib/ && cd htslib-1.3.2 && make && make install && cd .. && apt-get install -y openjdk-8-jre && echo "export PATH=/htslib-1.3.2:${PATH}" >> /root/.bashrc && echo "export LD_LIBRARY_PATH=/htslib-1.3.2:${PATH}" >> /root/.bashrc && bash install_2.sh && apt-get install -y libz-dev && apt-get install -y gcc  && apt-get install -y libgsl-dev && apt-get install -y bc && bash install_packages_2.sh  && python3 manage.py makemigrations && python3 manage.py migrate && Rscript install_tidyverse.R && wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2 && tar -xjvf bcftools.tar.bz2 && cd bcftools-1.3.1 && make && sudo make prefix=/usr/local/bin install && sudo ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools && cd .. && pip3 install tagore && sudo apt-get install -y librsvg2-bin && mkdir /usr/lib/tagore-data/ &&cp base.svg.p /usr/lib/tagore-data/base.svg.p && sudo apt-get remove cmake && sudo apt-get install -y build-essential && wget https://github.com/Kitware/CMake/releases/download/v3.20.5/cmake-3.20.5.tar.gz && tar xf cmake-3.20.5.tar.gz && cd cmake-3.20.5 && sudo apt-get install -y libssl-dev && ./configure && make && sudo apt-get install -y checkinstall && sudo checkinstall -y && cd .. && cd XGMix-master/ && pip3 install -r requirements.txt && cd ..


RUN apt-get update && apt-get install -y python3.7 && apt-get install -y python3-pip && pip3 install -r requirements.txt
#&& && apt-get install -y make  &&  apt-get install -y r-base && python3 manage.py makemigrations && python3 manage.py migrate

#RUN apt-get update && apt-get install -y  python3.7 && apt-get install -y python3-pip && apt autoremove && apt autoclean && apt-get clean packages && pip3 install -r requirements.txt && apt autoremove && apt autoclean && apt-get clean packages && apt-get install -y curl && apt-get install -y wget && wget https://github.com/wkhtmltopdf/wkhtmltopdf/releases/download/0.12.4/wkhtmltox-0.12.4_linux-generic-amd64.tar.xz && tar -xvf wkhtmltox-0.12.4_linux-generic-amd64.tar.xz && mv wkhtmltox/bin/wkhtmltopdf /usr/bin/wkhtmltopdf &&  mv wkhtmltox/bin/wkhtmltoimage /usr/bin/wkhtmltoimage && apt-get install -y libxrender1 && apt-get install -y libfontconfig && apt-get install -y libxext6 && apt-get install -y libssl1.0-dev && DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && apt-get install -y r-base && apt-get install -y sudo && bash install.sh && apt-get install -y make && apt-get install -y libz-dev && apt-get install -y gcc && apt-get install -y libgsl-dev && apt-get install -y bc && tar -xjvf htslib.tar.bz2 && export LD_LIBRARY_PATH=/usr/local/lib/ && cd htslib-1.3.2 && make && make install && cd .. && apt-get install -y openjdk-8-jre && echo "export PATH=/htslib-1.3.2:${PATH}" >> /root/.bashrc && echo "export LD_LIBRARY_PATH=/htslib-1.3.2:${PATH}" >> /root/.bashrc &&  apt-get install -y r-base && bash install_packages_2.sh && python3 manage.py makemigrations && python3 manage.py migrate
#RUN apt-get update && apt-get install -y  python3.7 && apt-get install -y python3-pip && apt autoremove && apt autoclean && apt-get clean packages && pip3 install -r requirements.txt && apt autoremove && apt autoclean && apt-get clean packages && apt-get install -y curl && apt-get install -y wget && wget https://github.com/wkhtmltopdf/wkhtmltopdf/releases/download/0.12.4/wkhtmltox-0.12.4_linux-generic-amd64.tar.xz && tar -xvf wkhtmltox-0.12.4_linux-generic-amd64.tar.xz && mv wkhtmltox/bin/wkhtmltopdf /usr/bin/wkhtmltopdf &&  mv wkhtmltox/bin/wkhtmltoimage /usr/bin/wkhtmltoimage && apt-get install -y libxrender1 && apt-get install -y libfontconfig && apt-get install -y libxext6 && apt-get install -y libssl1.0-dev && python3 manage.py makemigrations && python3 manage.py migrate && python3 manage.py runserver
