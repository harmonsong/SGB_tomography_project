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


def submit_job(source_this,jobs_status,dir_origin,email_content):
    # submit job
    try:
        os.chdir(dir_origin+'/src'+str(source_this)+'/')
        shell_order = "bsub -J 'src" + str(source_this) +  "' < runALL.sh"
        results = subprocess.run(shell_order, shell=True, check=True, capture_output=True,text=True)
        output = results.stdout.strip()
        job_id = output.split()[1][1:-1]
        print("Job submitted successfully! Job ID: ", str(job_id), '; source: ', str(source_this))
        email_content.append('Job submitted successfully! Job ID: '+str(job_id)+ '; source: '+ str(source_this))
        jobs_status[source_this] = 1
    except subprocess.CalledProcessError as e:
        email_content.append('Job finished fail! Job ID: '+str(job_id)+ ' source: '+ str(source_this))
        print("Job finished fail! Job ID: ", str(job_id), ' source: ', str(source_this))
        jobs_status[source_this] = -1
    return job_id, jobs_status, email_content

def check_status(jobs_status):
    global flag_max
    global flag_jobs
    global source_this
    global source_end
    global jobs_ids
    global dir_origin

    email_content = []
    # check if status has changed
    for job in jobs_ids.keys():
        job_id = jobs_ids[job]
        if jobs_status[job] < 0 or jobs_status[job] > 2:
            continue
        results = subprocess.run('bjobs '+job_id, shell=True, check=True, capture_output=True,text=True)
        output = results.stdout.strip()
        if 'not' in output:
            jobs_status[job] = 3
            email_content.append('Job finished successfully! Job ID: '+str(job_id), + ' source: '+ str(job))
            print("Job finished successfully! Job ID: ", str(job_id),' source: ', str(job))
        else:
            state = output.split()[10]

        if state == 'DONE':
            jobs_status[job] = 3
            email_content.append('Job finished successfully! Job ID: '+str(job_id)+ ' source: '+ str(job))
            print("Job finished successfully! Job ID: ", str(job_id), ' source: ', str(job))
        elif state == 'EXIT':
            jobs_status[job] = -1
            email_content.append('Job finished fail! Job ID: '+str(job_id) + ' source: '+ str(job))
            print("Job finished fail! Job ID: ", str(job_id), ' source: ', str(job))
        elif state == 'RUN' and jobs_status[job] != 2:
            jobs_status[job] = 2
            email_content.append('Job start running! Job ID: '+str(job_id) + ' source: '+ str(job))
            print("Job start running! Job ID: ", str(job_id), ' source: ', str(job))
    # submit jobs if there are spaces
    if email_content != []:
        time.sleep(5)
    status  = np.array(list(jobs_status.values()))
    flag_pend = np.where(status==1)[0]
    while len(flag_pend) < flag_max and source_this <= source_end:
        job_id,jobs_status,email_content = submit_job(source_this,jobs_status,dir_origin,email_content)
        jobs_ids[source_this] = job_id
        source_this += 1
        status  = np.array(list(jobs_status.values()))
        flag_pend = np.where(status==1)[0]
        time.sleep(5)
    return email_content,jobs_status
        


# parameters
flag_email = 1
flag_max = 10 # maximum pending jobs
source_start = 64 # starting source num
source_end = 66 # ending source num
jobs_status = {} 
jobs_ids = {}

# initialization
    # add previous pending jobs
results = subprocess.run('bjobs ', shell=True, check=True, capture_output=True,text=True)
output = results.stdout.strip().split()
pend_id_index = np.where(np.array(output) == 'PEND')[0]-2
pend_name_index = np.where(np.array(output) == 'PEND')[0]+3
pend_ids = np.array(output)[list(pend_id_index)]
pend_names = np.array(output)[list(pend_name_index)]
for i in range(len(pend_ids)):
    job_name = pend_names[i]
    jobs_ids[job_name] = pend_ids[i]
    jobs_status[job_name] = 1
    print('Previous pending job exist! Job ID: ', str(pend_ids[i]), '; job name: ', job_name)
    # add current aiming jobs
for source in range(source_start,source_end+1):
    jobs_status[source_start] = 0
source_this = source_start

# submit jobs
while 0 in list(jobs_status.values()) or 1 in list(jobs_status.values()) or 2 in list(jobs_status.values()):
    email_content,jobs_status = check_status(jobs_status)
    if email_content != [] and flag_email:
        email_to = [user_recive]
        email_title = 'Taiyi RaWave Jobs status change'
        email_content = '\n'.join(email_content)
        yag_server.send(email_to, email_title, email_content)
    time.sleep(20)
yag_server.close()


