#import urllib2,json
import sys
import os
import subprocess
import os.path

input_file = sys.argv[1]
multiple_files = "true"

path_prefix = "/scripts/"

print(len(sys.argv))
if(len(sys.argv) > 2):
    multiple_files_str = sys.argv[2]
    print(multiple_files_str)
    if(multiple_files_str == "multiple_files"):
        multiple_files = "true"


if(multiple_files == "true"):
    files_str = input_file
    files_list = files_str.split("_SEPERATOR_")
    for file_curr in files_list:
        #file_curr = path_prefix + file_curr_2
        print(file_curr)
        if not(os.path.isfile(file_curr)):
            continue
        file_1 = open(file_curr)
        htmlstr=file_1.read()
        file_1.close()
        curl_post = "curl --header \"Content-Type: application/json\" --request POST --data '" + str(htmlstr) + "' https://us-central1-enigmagenomics.cloudfunctions.net/raw_reports?action=create"
        output = subprocess.check_output(curl_post, shell=True)
        print(output)
else:        
    #input_file_new = path_prefix + input_file
    if not(os.path.isfile(input_file)):
        os._exit(1)
        #sys.exit("not successful - file not found")
    file_1 = open(input_file)
    htmlstr=file_1.read()
    file_1.close()
    curl_post = "curl --header \"Content-Type: application/json\" --request POST --data '" + str(htmlstr) + "' https://us-central1-enigmagenomics.cloudfunctions.net/raw_reports?action=create"
    output = subprocess.check_output(curl_post, shell=True)
    print(output)
