import subprocess
import yagmail
import sys
import time
import numpy as np
import os

# email config
user_send = '2947438278@qq.com'
user_recive = '2947438278@qq.com'
password = 'wmmnxrzfusqbddbi'
smtp_server = 'smtp.qq.com'
yag_server = yagmail.SMTP(user= user_send, password= password, host= smtp_server)
dir_origin = os.getcwd()


def submit_job(job_status,dir_origin,email_content):
    # submit job
    try:
        #os.chdir(dir_origin+'/template/')
        shell_order = "bsub -J 'template' < runALL.sh"
        results = subprocess.run(shell_order, shell=True, check=True, capture_output=True,text=True)
        output = results.stdout.strip()
        job_id = output.split()[1][1:-1]
        print("Job submitted successfully! Job ID: ", str(job_id), '; template')
        email_content.append('Job submitted successfully! Job ID: '+str(job_id)+ '; template')
        job_status = 1
    except subprocess.CalledProcessError as e:
        email_content.append('Job finished fail! Job ID: '+str(job_id)+'; template')
        print("Job finished fail! Job ID: ", str(job_id), ' ; template')
        job_status = -1
    return job_id, job_status, email_content

def check_status(job_status):
    global job_id
    global dir_origin

    email_content = []
    # submit jobs if there are spaces
    if job_status == 0:
        job_id,job_status,email_content = submit_job(job_status,dir_origin,email_content)
    # check if status has changed

    results = subprocess.run('bjobs '+job_id, shell=True, check=True, capture_output=True,text=True)
    output = results.stdout.strip()
    if 'not' in output:
        job_status = 3
        email_content.append('Job finished successfully! Job ID: '+str(job_id), + '; template')
        print("Job finished successfully! Job ID: ", str(job_id),'; template')
    else:
        state = output.split()[10]

    if state == 'DONE':
        job_status = 3
        email_content.append('Job finished successfully! Job ID: '+str(job_id)+ '; template')
        print("Job finished successfully! Job ID: ", str(job_id), '; template')
    elif state == 'EXIT':
        job_status = -1
        email_content.append('Job finished fail! Job ID: '+str(job_id) + '; template')
        print("Job finished fail! Job ID: ", str(job_id), '; template')
    elif state == 'RUN' and job_status != 2:
        job_status = 2
        email_content.append('Job start running! Job ID: '+str(job_id) + '; template')
        print("Job start running! Job ID: ", str(job_id),'; template')
    return email_content,job_status


flag_email = 1
job_status = 0 
job_id = 0

while job_status == 0 or job_status == 1 or job_status ==2:
    email_content,job_status = check_status(job_status)
    if email_content != [] and flag_email:
        email_to = [user_recive]
        email_title = 'Taiyi RaWave Jobs status change'
        email_content = '\n'.join(email_content)
        yag_server.send(email_to, email_title, email_content)
    time.sleep(60)
yag_server.close()