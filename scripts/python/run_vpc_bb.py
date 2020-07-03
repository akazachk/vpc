####
## First argument: the path to an mps/lp file or to an instances/batch file ([file].batch or [file].instances)
## Second argument: full path to output directory (accepted and necessary if the first argument is an instance, not batch/instances file; for a batch/instances file, we will specify the output directory stub in the first line; see below) 
##
## The first line of a .batch or .instances file will be the directory ``stub'',
## and output will be sent to ${PROJ_DIR}/results/stub if batch mode is off,
## and to ${PROJ_DIR}/results/stub/batchname if it is on.
## Each line contains either a relative input path to an instance (without the extension) or a batch name.
## The path is relative to ${PROJ_DIR}/data/instances, e.g., the line will be original/miplib2/bm23.
## A batch name is distinguished by having the batch end with a '/', e.g., '2/' or 'batch2/'.

## Set up proper path variables
import os
PROJ_DIR = os.path.abspath(os.environ['VPC_DIR'])
EXECUTABLE = PROJ_DIR + "/Release/vpc"
CUT_TYPE = 'vpc'
instances_path = os.path.abspath(os.environ['INSTANCE_DIR'])

## Solver options
depthList = [2,4,8,16,32,64]

## Set up output and input folders
results_path = PROJ_DIR + '/results'
outinfo_stub = CUT_TYPE + '-bb' 

## Get arguments
from sys import argv
assert (len(argv) > 1)
inst = os.path.abspath(argv[1])
is_instance = (inst[-4:] == '.mps') or (inst[-3:] == '.lp') or (inst[-7:] == '.mps.gz') or (inst[-6:] == '.lp.gz') or (inst[-8:] == '.mps.bz2') or (inst[-7:] == '.lp.bz2')
if is_instance:
    if (len(argv) > 1):
        results_path = argv[2] 
    else:
        raise Exception('If instance is given, so too should the full path to the results directory.')

## Where are the instances?
if is_instance:
    list_to_use = [inst]
else:
    with open(inst) as f_in:
      list_to_use = list(filter(None, (line.rstrip() for line in f_in)))

    ## The first line will be the name of the directory we should use
    results_path = results_path + '/' + list_to_use[0]
    list_to_use = list_to_use[1:]

## Finalize outinfo
#outinfo_dir = results_path + ('/' + dir_stub if len(dir_stuB) > 0 else "")
os.system("mkdir -p " + results_path)  # make the dir if it does not exist

## Choose order so that deepest for loop are the results you want to see first, fixing all others
batch_name = ''
for depth in depthList:
  for inst in list_to_use:
    ## Check if batch name
    if (inst[-1] == '/'):
      batch_name = inst
      continue

    ## Check if need to add "mps"
    inst_name = inst
    if (inst[-4:] != '.mps') and (inst[-3:] != '.lp') and (inst[-7:] != '.mps.gz') and (inst[-6:] != '.lp.gz') and (inst[-8:] != '.mps.bz2') and (inst[-7:] != '.lp.bz2'):
      inst_name = inst_name + '.mps'

    ## Run on instances_path/inst.mps
    infile = instances_path + '/' + inst_name
    curr_out_dir = results_path + '/' + batch_name
    outinfo = curr_out_dir + outinfo_stub + ".csv"
    #inst_stub = inst_name.split('/')[-1]
    #inst_stub = inst_stub.split('.')[0]  # problems if inst name has periods
    #paramfile = paramfile_dir + "/" + inst_stub + "_params.txt"
    
    ## In case the out directory does not exist
    os.system("mkdir -p " + curr_out_dir)

    ## Arguments
    extraparams = ' --optfile=' + PROJ_DIR + '/data/ip_obj.csv'
    extraparams = extraparams + ' --rounds=1'
    extraparams = extraparams + ' -d ' + str(depth)
    extraparams = extraparams + ' -t 3600'
    extraparams = extraparams + ' --bb_runs=1' 
    extraparams = extraparams + ' --bb_mode=10'
    extraparams = extraparams + ' --use_all_ones=1'
    extraparams = extraparams + ' --use_iter_bilinear=1'
    extraparams = extraparams + ' --use_disj_lb=1'
    extraparams = extraparams + ' --use_tight_points=0'
    extraparams = extraparams + ' --use_tight_rays=0'
    extraparams = extraparams + ' --use_unit_vectors=0'
    extraparams = extraparams + ' --gomory=-1'
    cmd = EXECUTABLE + ' -f ' + infile + ' --logfile=' + outinfo + extraparams 
    print(cmd)
    os.system(cmd + " > /dev/null 2>&1") 
