
This docker container includes the scripts for disease gene extraction, polygenic risk score calculation, global ancestry and local ancestry inference.

Build the docker container by running 

sudo docker-compose up

in this directory.

After doing that, you can log in as Admin (username: admin, pw: 1234) or as User (username: testuser, pw:1234) and use the following pages:

dis_calc/login.html - log in page.
dis_calc/overview_for_uploads.html - upload form for admin.
dis_calc/overview_for_customers.html - overview of result pages for customers.

You can find test files in the folder testfiles. The file for PRS can also be used for ancestry inference.

The results are written to json files that are linked in the result pages. The local ancestry results are stored as a png file.
