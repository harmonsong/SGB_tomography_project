import subprocess
import yagmail
import sys
import time
import numpy as np
import os

# email config
user_send = '2947438278@qq.com'
user_recive = '2947438278@qq.com'
password = 'mzltobvihgeidgih'
smtp_server = 'smtp.qq.com'
yag_server = yagmail.SMTP(user= user_send, password= password, host= smtp_server)
dir_origin = os.getcwd()


def submit_job(user_recive, yag_server,flag_this,jobs_status,dir_origin):
    # submit job
    try:
        os.chdir(dir_origin+'/src'+str(flag_this+1)+'/')
        shell_order = 'bsub < runALL.sh'
        results = subprocess.run(shell_order, shell=True, check=True, capture_output=True,text=True)
        output = results.stdout.strip()
        job_id = output.split()[1][1:-1]
        print("Job submitted successfully! Job ID: ", job_id, ' source: ', str(flag_this+1))
        email_to = [user_recive]
        email_title = 'Taiyi Job'
        email_content = 'Job submitted successfully! Job ID: '+job_id+ ' source: '+ str(flag_this+1)
        yag_server.send(email_to, email_title, email_content)
        jobs_status[flag_this] = 1
    except subprocess.CalledProcessError as e:
        print("Job submitted fail!: ", e)
        email_to = [user_recive]
        email_title = 'TaiYi Job'
        email_content = 'Job finished fail! Job ID: '+job_id+ ' source: '+ str(flag_this+1)
        yag_server.send(email_to, email_title, email_content)
        print("Job finished fail! Job ID: ", job_id, ' source: ', str(flag_this+1))
        jobs_status[flag_this] = -1
    return job_id, jobs_status

def check_status(jobs_status,jobs_ids):
    # enumerate running jobs
    flag_runs = np.where(jobs_status==1)[0]
    for flag_run in flag_runs:
        job_id = jobs_ids[flag_run]
        results = subprocess.run('bjobs '+job_id, shell=True, check=True, capture_output=True,text=True)
        output = results.stdout.strip()
        if 'not found' in output:
            email_to = [user_recive]
            email_title = 'TaiYi Job'
            email_content = 'Job finished successfully! Job ID: '+job_id + ' source: '+ str(flag_run+1)
            yag_server.send(email_to, email_title, email_content)
            print("Job finished successfully! Job ID: ", job_id, ' source: ', str(flag_run+1))
            jobs_status[flag_run] = 3
        else:
            state = output.split()[10]
            if state == 'DONE':
                email_to = [user_recive]
                email_title = 'TaiYi Job'
                email_content = 'Job finished successfully! Job ID: '+job_id + ' source: '+ str(flag_run+1)
                yag_server.send(email_to, email_title, email_content)
                print("Job finished successfully! Job ID: ", job_id, ' source: ', str(flag_run1))
                jobs_status[flag_run] = 3
            elif state == 'EXIT':
                email_to = [user_recive]
                email_title = 'TaiYi Job'
                email_content = 'Job finished fail! Job ID: '+job_id + ' source: '+ str(flag_run)
                yag_server.send(email_to, email_title, email_content)
                print("Job finished fail! Job ID: ", job_id, ' source: ', str(flag_run))
                jobs_status[flag_run] = -1
            elif state == 'RUN' and jobs_status[flag_run] == 1:
                email_to = [user_recive]
                email_title = 'TaiYi Job'
                email_content = 'Job Start running! Job ID: '+job_id
                yag_server.send(email_to, email_title, email_content)
                print("Job start running! Job ID: ", job_id, ' source: ', str(flag_run))
                jobs_status[flag_run] = 2
    return jobs_status

flag_max = 10
flag_jobs = 63 #63
jobs_status = np.zeros(flag_jobs) # 0--not start, 1--pend, 2--run, 3--finished, -1--fail
jobs_ids = {}

flag_this = 0
while flag_this < flag_jobs:
    # submit job
    job_id,jobs_status = submit_job(user_recive, yag_server,flag_this,jobs_status,dir_origin)
    jobs_ids[flag_this] = job_id
    if len(jobs_status[jobs_status==1]) < flag_max:
        flag_this += 1
        continue
    # check job status
    jobs_status = check_status(jobs_status,jobs_ids)
    while len(jobs_status[jobs_status==1]) == flag_max:
        time.sleep(600)
        jobs_status = check_status(jobs_status,jobs_ids)
    flag_this += 1

yag_server.close()