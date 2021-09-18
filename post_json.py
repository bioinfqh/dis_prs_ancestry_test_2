#import urllib2,json
import sys
import os
import subprocess

input_file = sys.argv[1]

file_1 = open(input_file)
htmlstr=file_1.read()
file_1.close()

curl_post = "curl --header \"Content-Type: application/json\" --request POST --data '" + str(htmlstr) + "' https://us-central1-enigmagenomics.cloudfunctions.net/raw_reports?action=create"


output = subprocess.check_output(curl_post, shell=True)

print(output)
