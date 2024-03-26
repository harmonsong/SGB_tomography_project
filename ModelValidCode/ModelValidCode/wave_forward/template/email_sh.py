import subprocess
import yagmail
import sys
import time

# shell name
shell_script_path = 'runALL.sh'

# email config
user_send = '2947438278@qq.com'
user_recive = '2947438278@qq.com'
password = 'mzltobvihgeidgih'
smtp_server = 'smtp.qq.com'
yag_server = yagmail.SMTP(user= user_send, password= password, host= smtp_server)

# submit job
try:
    results = subprocess.run('bsub < '+shell_script_path, shell=True, check=True, capture_output=True,text=True)
    output = results.stdout.strip()
    job_id = output.split()[1][1:-1]
    print("Job submitted successfully! Job ID: ", job_id)
    email_to = [user_recive]
    email_title = 'Taiyi Job'
    email_content = 'Job submitted successfully! Job ID: '+job_id
    yag_server.send(email_to, email_title, email_content)
except subprocess.CalledProcessError as e:
    print("Job submitted fail!: ", e)

# check job start 
flag = 0
while flag == 0:
    results = subprocess.run('bjobs '+job_id, shell=True, check=True, capture_output=True,text=True)
    output = results.stdout.strip()
    state = output.split()[10]
    if state == 'RUN':
        flag = 1
        email_to = [user_recive]
        email_title = 'TaiYi Job'
        email_content = 'Job Start running! Job ID: '+job_id
        yag_server.send(email_to, email_title, email_content)
        print("Job start running! Job ID: ", job_id)
    elif state == 'EXIT':
        flag = -1
        email_to = [user_recive]
        email_title = 'TaiYi Job'
        email_content = 'Job finished fail! Job ID: '+job_id
        yag_server.send(email_to, email_title, email_content)
        print("Job finished fail! Job ID: ", job_id)
        yag_server.close()
        sys.exit(1)
    else:
        time.sleep(120)


# check job finish
while flag == 1:
    results = subprocess.run('bjobs '+job_id, shell=True, check=True, capture_output=True,text=True)
    output = results.stdout.strip()
    if 'not found' in output:
        email_to = [user_recive]
        email_title = 'TaiYi Job'
        email_content = 'Job finished successfully! Job ID: '+job_id
        yag_server.send(email_to, email_title, email_content)
        print("Job finished successfully! Job ID: ", job_id)
        flag = 2
        break
    state = output.split()[10]
    if state == 'DONE':
        flag = 2
        email_to = [user_recive]
        email_title = 'Taiyi Job'
        email_content = 'Job finished successfully! Job ID: '+job_id
        yag_server.send(email_to, email_title, email_content)
        print("Job finished successfully! Job ID: ", job_id)
    elif state == 'EXIT':
        flag = -1
        email_to = [user_recive]
        email_title = 'Taiyi Job'
        email_content = 'Job finished fail! Job ID: '+job_id
        yag_server.send(email_to, email_title, email_content)
        print("Job finished fail! Job ID: ", job_id)
        yag_server.close()
        sys.exit(1)
    else:
        time.sleep(120)


yag_server.close()