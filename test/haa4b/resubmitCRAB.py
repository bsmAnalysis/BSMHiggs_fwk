import requests
import json
import os
import time
import subprocess
import getpass
import sys
import re

############################################## Variables and Functions Decleration Starts
runningMode = "" # 1 means whenever there's failed jobs, we will resubmit. If it's 0, we will only resubmit if this task has been finished, meaning only failed and finished jobs.
                # empty means we just check tasks, don't resubmit any tasks
#searchrange = 'from=&to=&timerange=lastWeek&pattern=' # Crab search range, you could change it to lastMonth
searchrange = 'from=2018-12-01+00%3A00&to=2019-01-11+23%3A59&timerange=&pattern=' # Crab search range, you could set it to your desired time range
def returnName(path,task):
  fileName = os.path.join(path,'crab.log')
  searchPattern = r'[0-9]{6}_[0-9]{6}:' + getpass.getuser() + '_' + task
  with open(fileName) as f:
    data =  ''.join(f.readlines())
  m = re.search(searchPattern, data)
  taskName = m.group()
  return taskName
############################################# Variables and Functions Decleration Stops 
  
if len(sys.argv) == 1:
  print("You must provide at least one input as the directory where you want to resubmit your crab jobs")
  print("If no second para, we will only check tasks. If second para is 1, we will resubmit any failed jobs; if it's 0, we resubmit tasks only with failed or finished tasks")
  print("For example: python resubmitCRAB.py results_2018_11_29 1")
  exit(1)

pathPre = os.getcwd()
path = os.path.join(os.path.join(pathPre,sys.argv[1]),"FARM/inputs")
if not os.path.isdir(path):
  print("Invalid path! Please make sure you input a directory and there's the folder '/FARM/inputs' under your input")
  exit(1)
DIRS = [name for name in os.listdir(path) if os.path.isdir(os.path.join(path,name)) and name.startswith('crab_')]
if len(DIRS)==0:
  print('Cannot find crab folder under the directory: {}. You cannot resubmit jobs if there is no crab folder!'.format(path))
  exit(1)

if len(sys.argv) > 2:
  runningMode = sys.argv[2]
  runningMode = int(runningMode)
  if runningMode == 1:
    print("We will resubmit tasks with failed jobs under the folder:\n {}".format(path))
  elif runningMode == 0:
    print("We will resubmit tasks only with failed or finished jobs under the folder:\n {}".format(path))
  else:
    print("Can not recognize the second parameter, only 0 or 1 allowed.")
    exit(1)
else:
  print("We will only check tasks under the folder: \n{}".format(path))

user = getpass.getuser()

url = 'http://dashb-cms-job.cern.ch/dashboard/request.py/antasktable?user={}&task=&{}'.format(user,searchrange)
headers = {
	'accept':'application/json, text/javascript, */*; q=0.01',
	'accept-Encoding':'gzip, deflate',
	'accept-Language':'zh-CN,zh;q=0.9,en;q=0.8',
	'connection':'keep-alive',
	'host':'dashb-cms-job.cern.ch',
#	'Referer':'http://dashb-cms-job.cern.ch/dashboard/templates/task-analysis/',
	'user-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.110 Safari/537.36',
#	'X-Requested-With':'XMLHttpRequest'
	'cache-Control':'max-age=0',
	'upgrade-Insecure-Requests':'1'
}
#params = {
#	'user':user,
#	'task':'',
#	'from':'',
#	'to':'',
#	'timerange':timerange,
#	'pattern':''
#}

#wbdata = requests.get(url,headers=headers,data=params,verify=False).text
wbdata = requests.get(url).text
data = json.loads(wbdata)

dirName = [] # a list to store jobs name
numFailed = [] # a list to store if there're any failed jobs
numFinished = [] # a list to store the number of finished jobs
numTotal = [] # a list to store the total jobs per task
flagFinished = [] # a list indicating if the task has been finished by CRAB, meanning there're only failed and finished jobs
flagCompleted = [] # a list indicating if the task has been completed by CRAB, meanning all jobs are finished
for job in data['antasks']:
  dirName.append(job['TASKNAME'])
  numFailed.append(job['FAILED'])
  numFinished.append(job['FINISHED'])
  numTotal.append(job['NJobTotal'])
  if (int(job['FAILED']) + int(job['FINISHED'])) == int(job['NJobTotal']):
    flagFinished.append(True)
  else:
    flagFinished.append(False)
  if int(job['FINISHED']) == int(job['NJobTotal']):
    flagCompleted.append(True)
  else:
    flagCompleted.append(False)

if runningMode == 1:
  flagFinished = [True] * len(flagFinished)
if runningMode == "":
  flagFinished = [False] * len(flagFinished)

count = 0
countFinished = 0 #counting how many jobs in the directory have been completed
totalTask = 0 # the number of tasks in total
totalFinished = 0 # the number of tasks finished
totalFound = 0 # the number of tasks we find on CRAB
errMsg = []
for task in DIRS:
  taskName = returnName(os.path.join(path,task),task)
  if taskName in dirName:
    index = dirName.index(taskName)
    totalFound += 1
    totalTask += numTotal[index]
    totalFinished += numFinished[index]
    if flagCompleted[index]:
      countFinished += 1
    else:
      print('Uncompleted Task: '+taskName)
    if flagFinished[index] and int(numFailed[index]) > 0: # to see if we have failed jobs for this task
#      print("{}. Resubmitting the job: {}, folder cteate: {}, task create: {}, delta: {}".format(count+1,task,ts,timeCreate[index],timeCreate[index]-ts))
      print("{}. Resubmitting the job: {}".format(count+1,task))
      command = "crab resubmit --dir=" + os.path.join(path,task)
#      subprocess.Popen(command,shell=True)
      result = subprocess.call(command,shell=True) # use call function to wait until the shell completes
#      if (not result == 0) and (count == 0):
      if not result == 0:
        print("Error occured when resubmitting jobs, please make sure you've sourced the crab script and get grid proxy certificates before running this script!")
        continue
#        exit(1)
      count += 1
  else:
    errMsg.append('Alert: can not find the task on CRAB: {}'.format(taskName))
if not totalFound == len(DIRS):
  for msg in errMsg:
    print(msg)
  print('-'*100)
  print('Alert: there are {} tasks in local but only {} tasks found on CRAB, see above for details'.format(len(DIRS),totalFound))
  print('If you submit jobs just now, you need to wait a few minutes to allow CRAB to process')
  print('Otherwise,if jobs are very old, you may want to enlarge the search range, say from lastWeek to lastMonth')
  print('Or maybe your submission failed for other reasons, please check with command "crab status"')
print('-'*100)
print('Summary Report')
print('Total {} tasks, finished: {}, resubmitted: {}'.format(len(DIRS),countFinished,count))
print('Total {} jobs, finished: {}, ratio: {}'.format(totalTask,totalFinished,float(totalFinished)/totalTask))
