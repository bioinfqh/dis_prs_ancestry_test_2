from django.shortcuts import render
from io import StringIO
from django.http import HttpResponse
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, redirect
from django.contrib.auth.hashers import check_password
from django.contrib.sites.shortcuts import get_current_site
from django.urls import reverse
from django.views import generic
from django.core.cache import cache
#from datetime import datetime
from django.contrib.auth.models import User
from django.shortcuts import render_to_response,render
from django.template import RequestContext
from django.contrib.auth import authenticate, login
from shutil import copyfile
import shutil
#### own imports
from django.contrib.auth import authenticate, login, logout
from dis_calc.make_report_from_files import make_prs
from dis_calc.make_report_from_files import make_dis_gene
from dis_calc.make_report_from_files import get_id_from_vcf
from dis_calc.make_report_from_files import make_ancestry_report
from dis_calc.make_report_from_files import make_local_anc_report
from django.core.cache import cache
import os.path
from django.core.mail import send_mail
from datetime import datetime

### *ACTUAL* imports (that have dependencies other than django and my own stuff) ####
import pandas as pd
import matplotlib.pyplot as plt

from dis_calc.models import user_with_data
from dis_calc.models import report
from dis_calc.models import ancestry_report
from dis_calc.models import local_ancestry_report
from django.contrib.auth.models import User, Group
import pdfkit 
import os


def logout_2(request):
	if(request.method == "POST"):
		logout(request)
		#return redirect('polls/logout.html')
		return redirect('dis_calc/logout.html')
	else:
		return render(request,'dis_calc/logout.html')

def login_2(request):
	if('username' in request.POST and 'password' in request.POST):
		username = request.POST['username']
		password = request.POST['password']
		user = authenticate(request, username=username, password=password)
		if user is not None:
			login(request, user)
			if user.groups.filter(name = "User").exists():
				return redirect('dis_calc/overview_for_customer.html')
			elif user.groups.filter(name = "Employee").exists():
				return redirect('dis_calc/overview_for_uploads.html')
				#return redirect('dis_calc/data_upload.html')
			#return render(request,'polls/login.html')
		else:
			text = "Username or password are incorrect"
			return render(request,'dis_calc/login.html',{'text':text})
			#return redirect('polls/clustering.html')
	else:
		return render(request,'dis_calc/login.html')
		#return redirect('polls/clustering.html')

# creates testuser and admin for testing.
def create_first_users(admin_name,customer_name):
    user1 = User.objects.create_user(admin_name, "bla@bla.com", "1234")
    user1.save()
    groups = Group.objects.all().values_list('name', flat=True)
    if not("Employee" in groups):
        group = Group(name = "Employee")
        group.save()
    else:
        group = Group.objects.get(name='Employee') 
    user1.groups.add(group)
    user2 = User.objects.create_user(customer_name, "bla@bla.com", "1234")
    user2.save()
    groups = Group.objects.all().values_list('name', flat=True)
    if not("User" in groups):
        group = Group(name = "User")
        group.save()
    else:
        group = Group.objects.get(name='User') 
    user2.groups.add(group)
    

    
def signup(request):
	if('username' in request.POST and 'password' in request.POST):
		username = request.POST['username']
		password = request.POST['password']
		email = ""
		if('email' in request.POST):
			email = request.POST['email']
		if(User.objects.filter(username=username).exists()):
			text = "Username already exists. Please choose another username!"
			return render(request,'dis_calc/signup.html',{'text':text})
		else:
			user = User.objects.create_user(username, email, password)
			user.save()
			if(username == "admin"):
				groups = Group.objects.all().values_list('name', flat=True)
				if not("Employee" in groups):
					group = Group(name = "Employee")
					group.save()
				else:
					group = Group.objects.get(name='Employee') 
				user.groups.add(group)
			else:
				groups = Group.objects.all().values_list('name', flat=True)
				if not("User" in groups):
					group = Group(name = "User")
					group.save()
				else:
					group = Group.objects.get(name='User') 
				user.groups.add(group)
			user_new = user_with_data(customer_id=username,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports="0",prs_report="",dis_report="")
			user_new.save()
			return render(request,'dis_calc/get_reports_customer.html',{'customer_id':username,'prs_report':"",'dis_gen_report':"",'multiple_reports':False})
			#return redirect('polls/clustering.html')
	else:
		return render(request,'dis_calc/signup.html')
		#return redirect('polls/clustering.html')


def make_pdf_from_str(htmlstr,outfile):
    fh=open("temp_out_file.html",'w+')
    fh.write(htmlstr)
    fh.close()
    pdfkit.from_file("temp_out_file.html", outfile)


# update PRS and disease gene report objects in database
def update_reports(username,report_types,disease_names,report_paths):
    reports = report.objects.filter(customer_id = username)
    # check if all parameter arrays are of the same length
    if not((len(report_types) == len(disease_names)) and (len(disease_names) == len(report_paths))):
        print("error")
        return("error")
    print("updating reports")
    for i in range(0,len(report_paths)):
        print(report_paths[i])
        report_type = report_types[i]
        disease_name = disease_names[i]
        # get all reports with is_current = yes
        report_curr = report.objects.filter(customer_id = username, is_current="yes",type_of_report=report_type,disease_type=disease_name)
        nbr_existing_reports = report.objects.filter(customer_id = username, is_current="yes",type_of_report=report_type,disease_type=disease_name).count()
        if(nbr_existing_reports == 0):
            continue
        for rep_curr in report_curr:
            print("setting is_current to no")
            # set is_current to no
            rep_curr.is_current="no"
            rep_curr.save()
        # generate new report object for new report
        report_new = report(customer_id=username,date=datetime.now(),is_current="yes",type_of_report=report_type,disease_type=disease_name,path=report_paths[i])
        report_new.save()
        if(report_type == "disease"):
            curr_user = user_with_data.objects.get(customer_id = username)
            # write names of analyzed diseases to user object in database
            dis_names_array = str(curr_user.diseases_with_dis_report).split(",")
            if not(disease_name in dis_names_array):
                dis_names_array.append(disease_name)
            dis_names_str = ",".join(dis_names_array)
            curr_user.diseases_with_dis_report=dis_names_str
            curr_user.save()


# update ancestry reports: set new report to is_current = yes and all other reports to is_current = no
def update_ancestry_reports(username,report_path,result_file_path):
    print("updating reports")
    #for rep_curr in reports:
    ## get all reports with is_current = yes
    report_curr = ancestry_report.objects.filter(customer_id = username, is_current="yes")
    nbr_existing_reports = ancestry_report.objects.filter(customer_id = username, is_current="yes").count()
    if not(nbr_existing_reports == 0):
        for rep_curr in report_curr:
            print("setting is_current to no")
            # set is_current to no
            rep_curr.is_current="no"
            rep_curr.save()
    # write new report object in database for new ancestry report
    report_new = ancestry_report(customer_id=username,date=datetime.now(),is_current="yes",result_path=result_file_path,path=report_path)
    report_new.save()
    



def update_local_ancestry_reports(username,img_path_new):
    print("updating reports")
    #for rep_curr in reports:
    ## get all reports with is_current = yes
    report_curr = local_ancestry_report.objects.filter(customer_id = username, is_current="yes")
    nbr_existing_reports = local_ancestry_report.objects.filter(customer_id = username, is_current="yes").count()
    if not(nbr_existing_reports == 0):
        for rep_curr in report_curr:
            print("setting is_current to no")
            # set is_current to no
            rep_curr.is_current="no"
            rep_curr.save()
    # write new report object in database for new ancestry report
    report_new = local_ancestry_report(customer_id=username,date=datetime.now(),is_current="yes",img_path=img_path_new)
    report_new.save()

def delete_user(request):
	if('username' in request.POST and 'password' in request.POST and request.user.is_authenticated):
		username = request.POST['username']
		password = request.POST['password']
		print(str(request.user))
		print(username)
		print(str(request.user.password))
		print(password)
		if(str(request.user) == username and check_password(password,request.user.password)):
			print("password found")
			u = User.objects.get(username=username)
			u.delete()
			user_dir = "user_uploaded_files/" + username
			if (os.path.isdir(user_dir)):
				shutil.rmtree(user_dir)
			text = "Your account was deleted."
			return render(request,'clustering/delete_user.html',{'text':text,'deleted':"true"})
		return render(request,'clustering/delete_user.html',{'text':""})
			#return redirect('polls/clustering.html')
	else:
		text = "Please input username and password!"
		return render(request,'clustering/delete_user.html',{'text':text})
		#return redirect('polls/clustering.html')



def errorpage(request):
		errors = ""
		if('errors' in request.session):
			errors = request.session['errors']
		errors_from_cache = cache.get('errors', 'has expired')
		if not(errors_from_cache == ""):
			errors = errors_from_cache
			cache.set('errors','')
		return render(request,'clustering/errorpage.html',{'errors':errors})


def pdf_test(request):
    if('input_file' in request.FILES):
        if(request.FILES['input_file']):
            html_str = request.FILES['input_file'].read().decode('utf-8','ignore')
            #fh=open("tmp_html.html",'w',encoding='utf-8', errors='ignore')
            #fh.write(html_str)
            #fh.close()
            make_pdf_from_str(html_str,"dis_calc/static/test_pdf.pdf")
        return render(request,'dis_calc/get_reports_customer.html',{'customer_id':"PDF_test",'multiple_reports':False,'prs_report':"test_pdf.pdf",'dis_gen_report':""})
    return render(request,'dis_calc/pdf_test.html')

def read_anc(path):
    table_file = open(path)
    htmlstr=table_file.read()
    table_file.close()
    lines = htmlstr.split("\n")
    anc_dict = {}
    for line in lines:
        if(len(line.split()) < 2):
            continue
        nat = line.split()[0]
        perc = line.split()[1]
        if("." in perc):
            anc_dict[nat] = float(perc)
    return(anc_dict)
        
### page with results for customers    
def result_page_customer(request):
    if request.user.is_authenticated:
        customer_id_curr = str(request.user)
        # calculate number of existing disease gene and PRS reports
        nbr_existing_reports = report.objects.filter(customer_id = customer_id_curr, is_current="yes",type_of_report="disease").count()
        nbr_existing_prs_reports = report.objects.filter(customer_id = customer_id_curr, is_current="yes",type_of_report="prs").count()
        print("Nbr PRS reports: " + str(nbr_existing_prs_reports))
        existing_reports = report.objects.filter(customer_id = customer_id_curr)
        print(existing_reports.all().values())
        existing_dis_reports = report.objects.filter(customer_id = customer_id_curr, is_current="yes",type_of_report="disease")
        ## get path of most recent PRS report (indicated by is_current)
        if(nbr_existing_prs_reports > 0):
            existing_prs_reports = report.objects.filter(customer_id = customer_id_curr, is_current="yes",type_of_report="prs")
            print(existing_prs_reports.all().values())
            prs_report = existing_prs_reports.first()
            prs_json_path = prs_report.path
            #print(prs_report_path)
        else:
            prs_json_path = None
        path_and_names = []
        print(existing_dis_reports.all().values())
        # get most recent disease gene report
        if(nbr_existing_reports == 1):
            dis_report = existing_dis_reports.first()
            dis_json_path = dis_report.path.replace("dis_calc/static/","")
        elif(nbr_existing_reports == 0):
            dis_report_path = None
        # generate list of disease name - path dictionaries for display of multiple disease reports on one html page. BY DEFAULT THIS IS NOT RUN YET
        if(nbr_existing_reports > 1):
            for report_curr in existing_dis_reports:
                dict_new = {}
                dict_new['disease_name'] = report_curr.disease_type
                dict_new['path'] = report_curr.path.replace("dis_calc/static/","")
            #dis_names = str(curr_user.diseases_with_dis_report).split(",")
            #report_names = dis_gen_report_name.split(",")
            #for i in range(0,len(dis_names)):
            #    dict_new['disease_name'] = dis_names[i]
            #    dict_new['path'] = report_names[i]
            #    path_and_names.append(dict_new)
            ## multiple_reports is true if a list of disease gene reports should be displayed, with paths and disease names indicated by path_and_names
            return render(request,'dis_calc/get_reports_customer.html',{'customer_id':customer_id_curr,'multiple_reports':True,'prs_json':prs_json_path,'path_and_names':path_and_names})
        return render(request,'dis_calc/get_reports_customer.html',{'customer_id':customer_id_curr,'multiple_reports':False,'prs_json':prs_json_path,'dis_gen_json':dis_json_path})
    else:
        return render(request, 'dis_calc/login.html')

# display of ancestry for customer
def ancestry_results_customer(request):
    if request.user.is_authenticated:
        customer_id_curr = str(request.user)
        # check if ancestry report exists for user
        nbr_existing_reports = ancestry_report.objects.filter(customer_id = customer_id_curr, is_current="yes").count()
        if(nbr_existing_reports > 0):
            existing_ancestry_reports = ancestry_report.objects.filter(customer_id = customer_id_curr, is_current="yes")
            ancestry_report_first = existing_ancestry_reports.first()
            anc_json_file_path = ancestry_report_first.path
            ancestry_result_path = ancestry_report_first.result_path
            print(anc_json_file_path)
            anc_dict = read_anc(ancestry_result_path)
            countries_and_percentages = []
            # write array of dicts with nationalities and ancestry percentages
            for nat in anc_dict:
                dict_new = {}
                dict_new['nat'] = nat
                dict_new['perc'] = anc_dict[nat]
                countries_and_percentages.append(dict_new)
        else:
            anc_json_file_path = None
        return render(request,'dis_calc/ancestry_overview_customer.html',{'customer_id':customer_id_curr,'ancestry_json':anc_json_file_path,'countries_and_percentages':countries_and_percentages})
    else:
        return render(request, 'dis_calc/login.html')


    
def local_ancestry_results_customer(request):
    if request.user.is_authenticated:
        customer_id_curr = str(request.user)
        # check if ancestry report exists for user
        nbr_existing_reports = local_ancestry_report.objects.filter(customer_id = customer_id_curr, is_current="yes").count()
        if(nbr_existing_reports > 0):
            existing_ancestry_reports = local_ancestry_report.objects.filter(customer_id = customer_id_curr, is_current="yes")
            ancestry_report_first = existing_ancestry_reports.first()
            anc_img_path = ancestry_report_first.img_path
        else:
            anc_img_path = None
        return render(request,'dis_calc/local_ancestry_report_customer.html',{'customer_id':customer_id_curr,'img_path':anc_img_path})
    else:
        return render(request, 'dis_calc/login.html')


def result_page_customer_OLD(request):
    if request.user.is_authenticated:
        customer_id_curr = str(request.user)
        curr_user = user_with_data.objects.get(customer_id = customer_id_curr)
        prs_report_name = curr_user.prs_report
        dis_gen_report_name = curr_user.dis_report
        nbr_reports = curr_user.has_reports
        if(float(nbr_reports) > 1.5):
            dis_names = str(curr_user.diseases_with_dis_report).split(",")
            report_names = dis_gen_report_name.split(",")
            for i in range(0,len(dis_names)):
                dict_new['disease_name'] = dis_names[i]
                dict_new['path'] = report_names[i]
                path_and_names.append(dict_new)
            return render(request,'dis_calc/get_reports_customer.html',{'customer_id':customer_id_curr,'multiple_reports':False,'prs_report':prs_report_name,'path_and_names':path_and_names,'dis_gen_report':dis_gen_report_name})
        return render(request,'dis_calc/get_reports_customer.html',{'customer_id':customer_id_curr,'multiple_reports':False,'prs_report':prs_report_name,'dis_gen_report':dis_gen_report_name})
    else:
        return render(request, 'dis_calc/login.html')
    
## overview of all uploads for employee
def overview_for_uploads(request):
    if(request.user.is_authenticated):
        if (request.user.groups.filter(name = "Employee").exists()):
            return render(request, 'dis_calc/overview_for_uploads.html')
    return render(request, 'dis_calc/login.html')
        
## overview page for customers        
def overview_for_customer(request):
    if(request.user.is_authenticated):
        if (request.user.groups.filter(name = "User").exists()):
            return render(request, 'dis_calc/overview_for_customer.html')
    return render(request, 'dis_calc/login.html')
        
### disease gene and/or PRS report generation for employees
def dis_calc_OLD(request):
    if(request.user.is_authenticated):
        if not (request.user.groups.filter(name = "Employee").exists()):
            return render(request, 'dis_calc/no_permission.html')
    else:
        return render(request, 'dis_calc/login.html')
    input_valid = "false"
    if('customer_id' in request.POST):
        customer_id_curr = request.POST['customer_id']
    else:
        customer_id_curr = "none"
    ## check if user has uploaded data
    if('prs_file' in request.FILES):
        if(request.FILES['prs_file']):
            input_valid = "true"
    elif('dis_file' in request.FILES):
        if(request.FILES['dis_file']):
            input_valid = "true"
    print(request.FILES)
    print(request.POST)
    print("starting method")
    print(input_valid)
    ## if no data were uploaded, assume the employee wants to upload new data
    if not(input_valid == "true"):
        return render(request,'dis_calc/data_upload.html')
    print("has valid data")
    prs_report_name = None
    if('prs_file' in request.FILES):
        if(request.FILES['prs_file']):
            # write PRS vcf to temporary file
            prs_vcf_str = request.FILES['prs_file'].read().decode('utf-8','ignore')
            fh=open("prs_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(prs_vcf_str)
            fh.close()
            bg_dataset="dis_calc/merged2filtered"
            patient_id = get_id_from_vcf("prs_vcf.vcf")
            # run PRS report in make_report_from_files
            make_prs(bg_dataset,"dis_calc/disease_list_5.txt","prs_vcf.vcf",patient_id,customer_id_curr)
            prs_report_name = "prs_report_" + customer_id_curr + ".pdf"
            dis_gen_report_name = None
            # generate entry in database for report
            report_new = report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",type_of_report="prs",disease_type="all",path=prs_report_name)
            report_new.save()
            user_count = user_with_data.objects.filter(customer_id = customer_id_curr).count()
            if(user_count > 0):
                curr_user = user_with_data.objects.filter(customer_id = customer_id_curr)
                curr_user.update(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=1,prs_report=prs_report_name,dis_report="",diseases_with_dis_report=(""))
            else:
                user_new = user_with_data(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=1,prs_report=prs_report_name,dis_report="",diseases_with_dis_report=(""))
                user_new.save()
            # update report objects in database
            update_reports(customer_id_curr,["prs"],["all"],[prs_report_name])
    #print(prs_vcf_str)
    if('dis_file' in request.FILES):
        if(request.FILES['dis_file']):
            # write input VCF to temporary file
            dis_vcf_str = request.FILES['dis_file'].read().decode('utf-8','ignore')
            fh=open("dis_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(dis_vcf_str)
            fh.close()
            patient_id=str(request.FILES['dis_file'].name).replace(".","")
            # make disease gene report
            [report_paths,disease_names_1] = make_dis_gene("dis_vcf.vcf",patient_id,"all",customer_id_curr)
            dis_gen_report_name = "dis_report_" + customer_id_curr + "_all.pdf"
            report_names = []
            report_paths = [i.replace("dis_calc/static/","") for i in report_paths]
            if(len(report_paths) > 0):
                report_names = report_paths
            report_path_str =  (",".join(report_names)).replace("dis_calc/static/","")
            user_count = user_with_data.objects.filter(customer_id = customer_id_curr).count()
            if(user_count > 0):
                curr_user = user_with_data.objects.filter(customer_id = customer_id_curr)
                curr_user.update(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=str(len(report_paths)),dis_report=report_path_str,diseases_with_dis_report=(",".join(disease_names_1)))
            else:
                user_new = user_with_data(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=str(len(report_paths)),dis_report=report_path_str,diseases_with_dis_report=(",".join(disease_names_1)))
                user_new.save()
            path_and_names = []
            dict_new = {}
            # update report objects in database (set current report to is_current = yes all other reports to is_current = no)
            update_reports(customer_id_curr,["disease"],["all"],[dis_gen_report_name])
            #generate report objects in database for all reports
            for i in range(0,len(report_paths)):
                dict_new['disease_name'] = disease_names_1[i]
                dict_new['path'] = report_paths[i].replace("dis_calc/static/","")
                path_and_names.append(dict_new)
                report_new =   report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",type_of_report="disease",disease_type=disease_names_1[i],path=report_paths[i])
                report_new.save()
                # update report objects in database
                update_reports(customer_id_curr,["disease"],[disease_names_1[i]],[report_paths[i].replace("dis_calc/static/","")])
            if(len(report_paths) > 1):
                return render(request,'dis_calc/get_reports.html',{'customer_id':customer_id_curr,'prs_report':prs_report_name,'path_and_names':path_and_names,'multiple_reports':True,'dis_gen_report':dis_gen_report_name})
    return render(request,'dis_calc/get_reports.html',{'customer_id':customer_id_curr,'prs_report':prs_report_name,'multiple_reports':False,'dis_gen_report':dis_gen_report_name})

def dis_calc(request):
    if(request.user.is_authenticated):
        if not (request.user.groups.filter(name = "Employee").exists()):
            return render(request, 'dis_calc/no_permission.html')
    else:
        return render(request, 'dis_calc/login.html')
    input_valid = "false"
    if('customer_id' in request.POST):
        customer_id_curr = request.POST['customer_id']
    else:
        customer_id_curr = "none"
    ## check if user has uploaded data
    if('prs_file' in request.FILES):
        if(request.FILES['prs_file']):
            input_valid = "true"
    elif('dis_file' in request.FILES):
        if(request.FILES['dis_file']):
            input_valid = "true"
    print(request.FILES)
    print(request.POST)
    print("starting method")
    print(input_valid)
    ## if no data were uploaded, assume the employee wants to upload new data
    if not(input_valid == "true"):
        return render(request,'dis_calc/data_upload.html')
    print("has valid data")
    prs_report_name = None
    prs_json_path = None
    dis_gen_json_path = None
    if('prs_file' in request.FILES):
        if(request.FILES['prs_file']):
            # write PRS vcf to temporary file
            prs_vcf_str = request.FILES['prs_file'].read().decode('utf-8','ignore')
            fh=open("prs_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(prs_vcf_str)
            fh.close()
            bg_dataset="dis_calc/merged2filtered"
            patient_id = get_id_from_vcf("prs_vcf.vcf")
            # run PRS report in make_report_from_files
            make_prs(bg_dataset,"dis_calc/disease_list_5.txt","prs_vcf.vcf",patient_id,customer_id_curr)
            prs_report_name = "prs_report_" + customer_id_curr + ".pdf"
            prs_json_path = "prs_" + customer_id_curr + ".json"
            dis_gen_report_name = None
            # generate entry in database for report
            report_new = report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",type_of_report="prs",disease_type="all",path=prs_json_path)
            report_new.save()
            user_count = user_with_data.objects.filter(customer_id = customer_id_curr).count()
            if(user_count > 0):
                curr_user = user_with_data.objects.filter(customer_id = customer_id_curr)
                curr_user.update(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=1,prs_report=prs_json_path,dis_report="",diseases_with_dis_report=(""))
            else:
                user_new = user_with_data(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=1,prs_report=prs_json_path,dis_report="",diseases_with_dis_report=(""))
                user_new.save()
            # update report objects in database
            update_reports(customer_id_curr,["prs"],["all"],[prs_json_path])
    #print(prs_vcf_str)
    if('dis_file' in request.FILES):
        if(request.FILES['dis_file']):
            # write input VCF to temporary file
            dis_vcf_str = request.FILES['dis_file'].read().decode('utf-8','ignore')
            fh=open("dis_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(dis_vcf_str)
            fh.close()
            patient_id=str(request.FILES['dis_file'].name).replace(".","")
            # make disease gene report
            [report_paths,disease_names_1] = make_dis_gene("dis_vcf.vcf",patient_id,"all",customer_id_curr)
            dis_gen_report_name = "dis_report_" + customer_id_curr + "_all.pdf"
            dis_gen_json_path = "dis_genes_" + customer_id_curr + "_all.json"
            report_names = []
            report_paths = [i.replace("dis_calc/static/","") for i in report_paths]
            if(len(report_paths) > 0):
                report_names = report_paths
            report_path_str =  (",".join(report_names)).replace("dis_calc/static/","")
            user_count = user_with_data.objects.filter(customer_id = customer_id_curr).count()
            if(user_count > 0):
                curr_user = user_with_data.objects.filter(customer_id = customer_id_curr)
                curr_user.update(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=str(len(report_paths)),dis_report=report_path_str,diseases_with_dis_report=(",".join(disease_names_1)))
            else:
                user_new = user_with_data(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=str(len(report_paths)),dis_report=report_path_str,diseases_with_dis_report=(",".join(disease_names_1)))
                user_new.save()
            path_and_names = []
            dict_new = {}
            # update report objects in database (set current report to is_current = yes all other reports to is_current = no)
            update_reports(customer_id_curr,["disease"],["all"],[dis_gen_json_path])
            #generate report objects in database for all reports
            for i in range(0,len(report_paths)):
                dict_new['disease_name'] = disease_names_1[i]
                dict_new['path'] = report_paths[i].replace("dis_calc/static/","")
                path_and_names.append(dict_new)
                report_new =   report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",type_of_report="disease",disease_type=disease_names_1[i],path=report_paths[i])
                report_new.save()
                # update report objects in database
                update_reports(customer_id_curr,["disease"],[disease_names_1[i]],[report_paths[i].replace("dis_calc/static/","")])
            if(len(report_paths) > 1):
                return render(request,'dis_calc/get_reports.html',{'customer_id':customer_id_curr,'prs_json':prs_json_path,'path_and_names':path_and_names,'multiple_reports':True,'dis_gen_report':dis_gen_report_name})
    return render(request,'dis_calc/get_reports.html',{'customer_id':customer_id_curr,'prs_json':prs_json_path,'multiple_reports':False,'dis_gen_json':dis_gen_json_path})


def ancestry_calc_OLD(request):
    if(request.user.is_authenticated):
        if not (request.user.groups.filter(name = "Employee").exists()):
            print("no permission")
            return render(request, 'dis_calc/no_permission.html')
    if not(request.user.is_authenticated):
        print("redirecting")
        return render(request, 'dis_calc/login.html')
    input_valid = "false"
    if('customer_id' in request.POST):
        customer_id_curr = request.POST['customer_id']
    else:
        customer_id_curr = "none"
    if('anc_file' in request.FILES):
        if(request.FILES['anc_file']):
            input_valid = "true"
    if not(input_valid == "true"):
        return render(request,'dis_calc/ancestry_upload.html')
    print("has valid data")
    if('anc_file' in request.FILES):
        if(request.FILES['anc_file']):
            # write input VCF content to temporary file
            anc_vcf_str = request.FILES['anc_file'].read().decode('utf-8','ignore')
            fh=open("dis_calc/anc_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(anc_vcf_str)
            fh.close()
            # make ancestry report
            [countries_and_percentages,outfile_for_ret,result_file] = make_ancestry_report("dis_calc/anc_vcf.vcf",customer_id_curr)
            # make report object in database
            report_new = ancestry_report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",result_path=result_file,path=outfile_for_ret)
            report_new.save()
            # update report objects in database
            update_ancestry_reports(customer_id_curr,outfile_for_ret,result_file)
            return render(request,'dis_calc/ancestry_overview.html',{'customer_id':customer_id_curr,'countries_and_percentages':countries_and_percentages,'ancestry_report_path':outfile_for_ret})
    return render(request,'dis_calc/ancestry_upload.html')
        
    
    
def ancestry_calc(request):
    if(request.user.is_authenticated):
        if not (request.user.groups.filter(name = "Employee").exists()):
            print("no permission")
            return render(request, 'dis_calc/no_permission.html')
    if not(request.user.is_authenticated):
        print("redirecting")
        return render(request, 'dis_calc/login.html')
    input_valid = "false"
    if('customer_id' in request.POST):
        customer_id_curr = request.POST['customer_id']
    else:
        customer_id_curr = "none"
    if('anc_file' in request.FILES):
        if(request.FILES['anc_file']):
            input_valid = "true"
    if not(input_valid == "true"):
        return render(request,'dis_calc/ancestry_upload.html')
    print("has valid data")
    if('anc_file' in request.FILES):
        if(request.FILES['anc_file']):
            # write input VCF content to temporary file
            anc_vcf_str = request.FILES['anc_file'].read().decode('utf-8','ignore')
            fh=open("dis_calc/anc_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(anc_vcf_str)
            fh.close()
            # make ancestry report
            [countries_and_percentages,outfile_for_ret,result_file] = make_ancestry_report("dis_calc/anc_vcf.vcf",customer_id_curr)
            anc_json_file_path = "ancestry_" + customer_id_curr + ".json"
            # make report object in database
            report_new = ancestry_report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",result_path=result_file,path=anc_json_file_path)
            report_new.save()
            # update report objects in database
            update_ancestry_reports(customer_id_curr,anc_json_file_path,result_file)
            return render(request,'dis_calc/ancestry_overview.html',{'customer_id':customer_id_curr,'countries_and_percentages':countries_and_percentages,'ancestry_report_path':outfile_for_ret,'ancestry_json':anc_json_file_path})
    return render(request,'dis_calc/ancestry_upload.html')
   
   
   
        
        
def local_ancestry_calc(request):
    if(request.user.is_authenticated):
        if not (request.user.groups.filter(name = "Employee").exists()):
            print("no permission")
            return render(request, 'dis_calc/no_permission.html')
    if not(request.user.is_authenticated):
        print("redirecting")
        return render(request, 'dis_calc/login.html')
    input_valid = "false"
    if('customer_id' in request.POST):
        customer_id_curr = request.POST['customer_id']
    else:
        customer_id_curr = "none"
    if('anc_file' in request.FILES):
        if(request.FILES['anc_file']):
            input_valid = "true"
    if not(input_valid == "true"):
        return render(request,'dis_calc/local_ancestry_upload.html')
    print("has valid data")
    if('anc_file' in request.FILES):
        if(request.FILES['anc_file']):
            # write input VCF content to temporary file
            anc_vcf_str = request.FILES['anc_file'].read().decode('utf-8','ignore')
            fh=open("dis_calc/local_anc_vcf.vcf",'w',encoding='utf-8', errors='ignore')
            fh.write(anc_vcf_str)
            fh.close()
            print("file written to local_anc_vcf.vcf")
            # make ancestry report
            anc_img = make_local_anc_report("dis_calc/local_anc_vcf.vcf",customer_id_curr)
            #[countries_and_percentages,outfile_for_ret,result_file] = make_ancestry_report("dis_calc/anc_vcf.vcf",customer_id_curr)
            # make report object in database
            report_new = local_ancestry_report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",img_path=anc_img)
            report_new.save()
            # update report objects in database
            update_local_ancestry_reports(customer_id_curr,anc_img)
            return render(request,'dis_calc/local_ancestry_report.html',{'customer_id':customer_id_curr,'img_path':anc_img})
    return render(request,'dis_calc/local_ancestry_upload.html')

    


def errorpage(request):
		errors = ""
		if('errors' in request.session):
			errors = request.session['errors']
		return render(request,'clustering/errorpage.html',{'errors':errors})




def infopage(request):			
	return render(request,'clustering/infopage.html')

def sources(request):			
	return render(request,'clustering/sources.html')



