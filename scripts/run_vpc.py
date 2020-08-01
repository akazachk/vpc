####
## Four arguments, first three are necessary
## First argument: the path to an mps/lp file or to an instances/batch file ([file].batch or [file].instances)
## Second argument: full path to output directory
## Third argument: running mode (bb, bb0, or preprocess)
## Fourth argument: optional extra args to pass to the solver
##
## The first line of a .batch file will be the directory "stub"
## and output will be sent to ${PROJ_DIR}/results/stub if batch mode is off,
## and to ${PROJ_DIR}/results/stub/batchname if it is on.
## Each line contains either a relative input path to an instance (without the extension) or a batch name.
## The path is relative to ${PROJ_DIR}/data/instances, e.g., the line will be original/miplib2/bm23.
## A batch name is distinguished by having the batch end with a '/', e.g., '2/' or 'batch2/'.

## Set up proper path variables
import os
from datetime import datetime
PROJ_DIR = os.path.abspath(os.environ['VPC_DIR'])
EXECUTABLE = PROJ_DIR + "/Release/vpc"
CUT_TYPE = 'vpc'
instances_path = os.path.abspath(os.environ['INSTANCE_DIR'])
results_path = PROJ_DIR + '/results'

## Script options
OPTION_TEST = -1
OPTION_BB = 0
OPTION_BB0 = 1
OPTION_PREPROCESS = 2
option = OPTION_BB

## Where are the instances?
from sys import argv
ARG_INST = 1
ARG_RESULT = 2
ARG_OPTION = 3
ARG_EXTRA_OPTS = 4
assert (len(argv) > 2)
try:
    option = argv[ARG_OPTION]
    if (option == "bb"):
        option = OPTION_BB
    elif option == "bb0":
        option = OPTION_BB0
    elif option == "preprocess":
        option = OPTION_PREPROCESS
    elif option == "test":
        option = OPTION_TEST
    else:
        option = int(option)
        assert (option == OPTION_BB or option == OPTION_BB0 or option == OPTION_PREPROCESS)
except:
    raise ValueError("Unable to process the first argument into one of the script options")
try:
    inst = os.path.abspath(argv[ARG_INST]) # make absolute path if given as relative
except:
    raise "Unable to read instances path from the second argument"
is_instance = (inst[-4:] == '.mps') or (inst[-3:] == '.lp') or (inst[-7:] == '.mps.gz') or (inst[-6:] == '.lp.gz') or (inst[-8:] == '.mps.bz2') or (inst[-7:] == '.lp.bz2')
try:
    results_path = argv[ARG_RESULT]
except:
    raise "Unable to read results path from the third argument"

## Solver options
depthList = [0]
outinfo_stub = CUT_TYPE + "-null"
extraparams = ' --optfile=' + PROJ_DIR + '/data/ip_obj.csv'
if option == OPTION_BB:
    depthList = [2,4,8,16,32,64]
    outinfo_stub = CUT_TYPE + '-bb'
    extraparams = extraparams + ' --rounds=1'
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
elif option == OPTION_BB0:
    outinfo_stub = CUT_TYPE + '-bb0'
    extraparams = extraparams + ' --rounds=1'
    extraparams = extraparams + ' -t 3600'
    extraparams = extraparams + ' --bb_runs=7'
    extraparams = extraparams + ' --bb_mode=001'
    extraparams = extraparams + ' --use_all_ones=1'
    extraparams = extraparams + ' --use_iter_bilinear=1'
    extraparams = extraparams + ' --use_disj_lb=1'
    extraparams = extraparams + ' --use_tight_points=0'
    extraparams = extraparams + ' --use_tight_rays=0'
    extraparams = extraparams + ' --use_unit_vectors=0'
elif option == OPTION_PREPROCESS:
    outinfo_stub = 'cleaning_log'
    extraparams = extraparams + ' --preprocess=1'
    extraparams = extraparams + ' -t 7200'
    extraparams = extraparams + ' --bb_runs=1'
    extraparams = extraparams + ' --bb_mode=001'
    # Set strategy (536 uses Gurobi, 532 uses CPLEX)
    extraparams = extraparams + ' --bb_strategy=532'
elif option == OPTION_TEST:
    depthList = [2]
    outinfo_stub = CUT_TYPE + '-test'
else:
    raise "Option not recognized"

userparams = ''
if (len(argv) >= 4):
    tmpuserparams = argv[ARG_EXTRA_OPTS]
    tmpuserparams = tmpuserparams.strip()
    if len(tmpuserparams) > 0:
        depthList = [0]
        userparams = ' ' + tmpuserparams

## Where to get instances
if is_instance:
    list_to_use = [inst]
else:
    with open(inst) as f_in:
      list_to_use = list(filter(None, (line.rstrip() for line in f_in)))

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
    infile = (instances_path + '/' if not is_instance else "") + inst_name
    curr_out_dir = results_path + '/' + batch_name
    outinfo = curr_out_dir + outinfo_stub + ".csv"
    #inst_stub = inst_name.split('/')[-1]
    #inst_stub = inst_stub.split('.')[0]  # problems if inst name has periods
    #paramfile = paramfile_dir + "/" + inst_stub + "_params.txt"

    ## In case the out directory does not exist
    os.system("mkdir -p " + curr_out_dir)

    ## Arguments
    cmd = EXECUTABLE + ' -f ' + infile + ' --logfile=' + outinfo + extraparams + ' -d ' + str(depth) + userparams
    now = datetime.now()
    print(now.strftime("%Y-%m-%d %H:%M:%S") + ": " + cmd, flush=True) # flush=True requires python 3.3+
    os.system(cmd + " > /dev/null 2>&1")
