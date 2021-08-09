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
from django.core.cache import cache
import os.path
from django.core.mail import send_mail
from datetime import datetime

### *ACTUAL* imports (that have dependencies other than django and my own stuff) ####
import pandas as pd
import matplotlib.pyplot as plt

from dis_calc.models import user_with_data
from dis_calc.models import report
from django.contrib.auth.models import User, Group

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
				return redirect('dis_calc/get_reports_customer.html')
			elif user.groups.filter(name = "Employee").exists():
				return redirect('dis_calc/data_upload.html')
			#return render(request,'polls/login.html')
		else:
			text = "Username or password are incorrect"
			return render(request,'dis_calc/login.html',{'text':text})
			#return redirect('polls/clustering.html')
	else:
		return render(request,'dis_calc/login.html')
		#return redirect('polls/clustering.html')

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



def result_page_customer(request):
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

def dis_calc(request):
    if(request.user.is_authenticated):
        if not (request.user.groups.filter(name = "Employee").exists()):
            return render(request, 'dis_calc/no_permission.html')
    else:
        return render(request, 'dis_calc/login.html')
    input_valid = "false"
    #print(request.POST)
    if('customer_id' in request.POST):
        customer_id_curr = request.POST['customer_id']
    else:
        customer_id_curr = "none"
    if('prs_file' in request.FILES and 'dis_file' in request.FILES):
        if(request.FILES['prs_file'] and request.FILES['dis_file']):
            input_valid = "true"
    print(request.FILES)
    print(request.POST)
    print("starting method")
    if not(input_valid == "true"):
        return render(request,'dis_calc/data_upload.html')
    print("has valid data")
    prs_vcf_str = request.FILES['prs_file'].read().decode('utf-8','ignore')
    #print(prs_vcf_str)
    dis_vcf_str = request.FILES['dis_file'].read().decode('utf-8','ignore')
    fh=open("prs_vcf.vcf",'w',encoding='utf-8', errors='ignore')
    fh.write(prs_vcf_str)
    fh.close()
    fh=open("dis_vcf.vcf",'w',encoding='utf-8', errors='ignore')
    fh.write(dis_vcf_str)
    fh.close()
    bg_dataset="dis_calc/merged2filtered"
    patient_id=str(request.FILES['dis_file'].name).replace(".","")
    print(patient_id)
    #make_prs(bg_dataset,"dis_calc/disease_list.txt","prs_vcf.vcf",patient_id,customer_id_curr)
    #make_prs(bg_dataset,"all","prs_vcf.vcf",patient_id,customer_id_curr)
    make_prs(bg_dataset,"dis_calc/disease_list_5.txt","prs_vcf.vcf",patient_id,customer_id_curr)
    #make_dis_gene("dis_vcf.vcf",patient_id,"dis_calc/disease_groups.txt",customer_id_curr)
    [report_paths,disease_names_1] = make_dis_gene("dis_vcf.vcf",patient_id,"all",customer_id_curr)
    #make_dis_gene("dis_vcf.vcf",patient_id,"dis_calc_application/dis_calc/static/disease_groups.txt",customer_id):
    #prs_report_name = "prs_report_" + customer_id_curr + ".pdf"
    #dis_gen_report_name = "dis_report_" + customer_id_curr + ".pdf"
    #prs_report_name = "your_prs_report.pdf"
    #dis_gen_report_name = "your_report.pdf"
    prs_report_name = "prs_report_" + customer_id_curr + ".pdf"
    dis_gen_report_name = "dis_report_" + customer_id_curr + "_all.pdf"
    report_names = []
    if(len(report_paths) > 0):
        report_names = report_paths
    report_path_str =  ",".join(report_names)
    user_count = user_with_data.objects.filter(customer_id = customer_id_curr).count()
    if(user_count > 0):
        curr_user = user_with_data.objects.filter(customer_id = customer_id_curr)
        curr_user.update(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=str(len(report_paths)),prs_report=prs_report_name,dis_report=report_path_str,diseases_with_dis_report=(",".join(disease_names_1)))
    else:
        user_new = user_with_data(customer_id=customer_id_curr,firstname="Max",lastname="Mustermann",email="bla@bla.com",has_reports=str(len(report_paths)),prs_report=prs_report_name,dis_report=report_path_str,diseases_with_dis_report=(",".join(disease_names_1)))
        user_new.save()
    report_new = report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",type_of_report="prs",disease_type="all",path=prs_report_name)
    report_new.save()
    #prs_report_name = "your_prs_report.pdf"
    #dis_gen_report_name = "your_report.pdf"
    path_and_names = []
    dict_new = {}
    for i in range(0,len(report_paths)):
        dict_new['disease_name'] = disease_names_1[i]
        dict_new['path'] = report_paths[i]
        path_and_names.append(dict_new)
        report_new = report(customer_id=customer_id_curr,date=datetime.now(),is_current="yes",type_of_report="disease",disease_type=disease_names_1[i],path=report_paths[i])
        report_new.save()
    if(len(report_paths) > 1):
        return render(request,'dis_calc/get_reports.html',{'customer_id':customer_id_curr,'prs_report':prs_report_name,'path_and_names':path_and_names,'multiple_reports':True,'dis_gen_report':dis_gen_report_name})
    return render(request,'dis_calc/get_reports.html',{'customer_id':customer_id_curr,'prs_report':prs_report_name,'multiple_reports':False,'dis_gen_report':dis_gen_report_name})





def errorpage(request):
		errors = ""
		if('errors' in request.session):
			errors = request.session['errors']
		return render(request,'clustering/errorpage.html',{'errors':errors})




def infopage(request):			
	return render(request,'clustering/infopage.html')

def sources(request):			
	return render(request,'clustering/sources.html')



