import subprocess
import yagmail
import sys
import time
import numpy as np
import os

def submit_job(pbs_name,flag_this,d_start,d_len,jobs_status,email_content):
    # submit job
    try:
        arg = 'dstart='+str(d_start)+',dlen='+str(d_len)
        shell_order = 'qsub -v '+ arg + ' ' + pbs_name
        results = subprocess.run(shell_order, shell=True, check=True, capture_output=True,text=True)
        output = results.stdout.strip()
        job_id = output.split()[0]
        print("Job submitted successfully! Job ID: ", job_id, '; day: ', str(d_start))
        email_content.append( 'Job submitted successfully! Job ID: '+job_id + '; day: '+str(d_start) )
        jobs_status[flag_this] = 1
    except subprocess.CalledProcessError as e:
        email_content.append('Job submit fail! Job ID: '+job_id + '; day: '+str(d_start))
        print("Job submit fail! Job ID: ", job_id, '; day: ', str(d_start))
        jobs_status[flag_this] = -1
    return job_id, jobs_status,email_content

def check_status(jobs_status):
    global flag_max
    global d_start0
    global d_len
    global flag_this

    global jobs_ids
    global pbs_name
    global user_name
    global flag_jobs
    global user_recive
    global yag_server

    email_content = []
    # enumerate running jobs
    flag_runs = np.where(jobs_status==1)[0]
    # submit jobs if there is spaces
    while len(flag_runs) < flag_max and flag_this < flag_jobs:
        job_id,jobs_status,email_content = submit_job(pbs_name,flag_this,d_start0+flag_this,d_len,jobs_status,email_content)
        jobs_ids[flag_this] = job_id
        flag_this += 1
        flag_runs = np.where(jobs_status==1)[0]

    # check if job has finised
    results = subprocess.run('qstat -u '+user_name, shell=True, check=True, capture_output=True,text=True)
    output = results.stdout.strip()
    for flag_run in flag_runs:
        job_id = jobs_ids[flag_run]
        if not job_id in output:
            print("Job finished successfully! Job ID: ", job_id, ', day: ', str(128+flag_run))
            jobs_status[flag_run] = 3
            email_content.append('Job finished successfully! Job ID: '+job_id+ ', day: '+ str(128+flag_run))
        elif jobs_status[flag_run] == 1:
            results1 = subprocess.run('qstat '+job_id, shell=True, check=True, capture_output=True,text=True)
            output1 = results1.stdout.strip()
            state = output1.split()[18]
            if state == 'R':
                print("Job start running! Job ID: "+job_id+', day: '+ str(128+flag_run))
                email_content.append("Job start running! Job ID: "+job_id+ ', day: '+str(128+flag_run))
                jobs_status[flag_run] = 2
            elif state == 'Q':
                continue
            else:
                print("Job finished successfully! Job ID: ", job_id, ', day: ', str(128+flag_run))
                jobs_status[flag_run] = 3
                email_content.append('Job finished successfully! Job ID: '+job_id+ ', day: '+ str(128+flag_run))
    return email_content,jobs_status


# email config
flag_email = 0 # 0--do not send email; 1--send email
user_send = '2947438278@qq.com'
user_recive = '2947438278@qq.com'
password = 'wmmnxrzfusqbddbi'
smtp_server = 'smtp.qq.com'
yag_server = yagmail.SMTP(user= user_send, password= password, host= smtp_server)

flag_max = 4
d_start0 = 128
d_len = 1
user_name = 'songsh'
pbs_name = '1-2-2.pbs'
flag_jobs = 31 #63
jobs_status = np.zeros(flag_jobs) # 0--not start, 1--pend, 2--run, 3--finished
jobs_ids = {}

flag_this = 0 # source_num -1
while 0 in jobs_status:
    email_content,jobs_status = check_status(jobs_status)
    if email_content != [] and flag_email:
        email_to = [user_recive]
        email_title = 'F-J correlation Job'
        email_content = '\n'.join(email_content)
        yag_server.send(email_to, email_title, email_content)
    time.sleep(5)
yag_server.close()