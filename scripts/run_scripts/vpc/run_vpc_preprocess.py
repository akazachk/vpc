####
## The inputted file is either [filename].instances or [file].batch.
## (The extension is not important.)
## The first line of this file will be the directory ``stub'',
## and output will be sent to ${PROJ_DIR}/results/instances/stub if batch mode is off,
## and to ${PROJ_DIR}/results/instances/batches/stub/batchname if it is on.
## Each line contains either a relative input path to an instance (without the extension) or a batch name.
## The path is relative to ${PROJ_DIR}/data/instances, e.g., the line will be miplib2/bm23.
## A batch name is distinguished by having the batch end with a '/', e.g., '2/' or 'batch2/'.

## Set up proper path variables
import os
PROJ_DIR = os.path.abspath(os.environ['VPC_DIR'])
EXECUTABLE = PROJ_DIR + "/Release/vpc"
CUT_TYPE = 'vpc'
instances_path = os.path.abspath(os.environ['INSTANCE_DIR'])

## Set up output and input folders
results_path = PROJ_DIR + '/results'
#paramfile_dir = PROJ_DIR + '/data/params'
#instances_path = PROJ_DIR + '/data/instances'
instances_file = instances_path + '/' + "test.instances"

outinfo_stub = "cleaning_log"
outinfo_dir = results_path

## Get arguments
from sys import argv
use_batches = False  # set to true/false depending on if mps files are all in one folder or divided up into subfolders
if (len(argv) > 1):
  use_batches = True if argv[1] in ['true', 'True', '1', 't'] else False
  if (use_batches and len(argv) < 2):
    raise ValueError('When using batches, specifying the folder is required')

if (len(argv) > 2):
  instances_file = os.path.abspath(argv[2])

## Where are the instances?
with open(instances_file) as f_in:
  list_to_use = list(filter(None, (line.rstrip() for line in f_in)))

## The first line will be the name of the directory we should use
dir_stub = list_to_use[0]
list_to_use = list_to_use[1:]

if use_batches:
  dir_stub = "batches/" + dir_stub

## Finalize outinfo
outinfo_dir = outinfo_dir + '/' + dir_stub 
os.system("mkdir -p " + outinfo_dir)  # make the dir if it does not exist

## Choose order so that deepest for loop are the results you want to see first, fixing all others
batch_name = ''
for inst in list_to_use:
  ## Check if batch name
  if (inst[-1] == '/'):
    batch_name = inst
    continue

  ## Check if need to add "mps"
  inst_name = inst
  if (inst[-4:] != '.mps') and (inst[-3:] != '.lp') and (inst[-7:] != '.mps.gz') and (inst[-6:] != '.lp.gz'):
    inst_name = inst_name + '.mps'

  ## Run on instances_path/inst.mps
  infile = instances_path + '/' + inst_name
  curr_out_dir = outinfo_dir + '/' + batch_name
  outinfo = curr_out_dir + outinfo_stub + ".csv"

  ## Get instance filename stub to use instance-specific parameters
  #inst_stub = inst_name.split('/')[-1]
  #inst_stub = inst_stub.split('.')[0]  # problems if inst name has periods
  #paramfile = paramfile_dir + "/" + inst_stub + "_params.txt"
  
  ## In case the out directory does not exist
  os.system("mkdir -p " + curr_out_dir)

  ## Arguments
  params = ' -f ' + infile
  params = params + ' --logfile=' + outinfo
  params = params + ' --optfile=' + PROJ_DIR + '/data/ip_obj.csv'
  params = params + ' --temp=1'
  params = params + ' -t 7200'
  params = params + ' --bb_runs=7' 
  params = params + ' --bb_mode=001'
  cmd = EXECUTABLE + params
  print(cmd)
  os.system(cmd + " > /dev/null 2>&1") 
