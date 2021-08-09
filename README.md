README:

This docker container contains a Django webserver with basic functionality for generating PRS and disease gene reports from VCF files, and to estimate ancestry based on VCF files.
You can run all necessary steps to create and run a Docker container with the server by running run_docker.sh.

The website contains the following pages:

login.html - Log in as user/admin. Existing logins are: testuser (pw: 1234) and admin (pw: 1234).
logout.html - Log out.
signup.html - Sign up as new user.
overview_for_uploads.html - Accessible for admin. Here you can select which type of analysis you want to run.
data_upload.html - Here you can upload files for a user and generate reports. This page will forward you to an overview of the user's reports once they are generated.
ancestry_upload.html - Here you can upload a VCF file for ancestry estimation.
overview_for_customer.html - Accessible for users (customers). Here you can select which reports you want to view.
get_reports_customer.html - Here all current disease gene and PRS reports of an user are listed.
ancestry_overview_customer.html - Here you can see your estimated ancestry. This page contains a link to an interactive world map where your ancestry is vizualised.

Data input:

Disease genes + PRS:

It is possible to upload a file for the disease risk reports, or a file for disease gene reports, or both.
The input format for both files is vcf. The vcf files for PRS calculation must include a customer name in the header column (which is always required in a standard vcf file). The files for disease gene extraction should be annotated VEP output.


Ancestry:

You can upload a VCF file.


For tests you can use the files provided at testfiles.
In the current setting, disease risk reports can be generated for the following diseases:
Breast Cancer, Cervical Cancer, Colorectal Cancer, Coronary Artery Disease, Hypertension, Prostate Cancer


Users:

Two user logins are pre-installed:

testuser (pw: 1234) - a simple customer that can access his reports.
admin (pw: 1234) - the admin (employee) that can upload files to generate reports.


About this container:
In the folder dis_calc there are two lists that specify which genes are included in analyses for certain diseases.
The list with diseases for risk reports is in a file called "disease_list_5.txt".
The list with disease groups for disease gene reports is in a file called "disease_groups.txt".
In a file called "gene_disease_groups.csv" all genes associated to each of the disease groups is listed.
For some calculations, synonyms of diseases must be collected and removed from output disease lists. There is a list of pre-compiled synonyms at "syndict_temp.txt". If an entry in a VEP output file lists a disease that is not listed in this file, new synonyms are collected from an API. This is the task that takes the longest for most analyses.
For better runtimes, only 20 lines of the VEP input file are used.
Ancestry calculation takes a few minutes.
