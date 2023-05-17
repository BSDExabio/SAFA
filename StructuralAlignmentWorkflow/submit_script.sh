#!/bin/bash

#SBATCH -A ACCOUNT
#SBATCH -t 36:00:00
#SBATCH -N 64
#SBATCH -p batch
#SBATCH -J usalign
#SBATCH -o out.%J
#SBATCH -e err.%J

# for ANDES batch queue, nSchedulerCores + nClientCores should never be greater than 32; ideally, this sum would consume all cores on the primary node OR leave some resources available for a worker to be spun up on the primary/head node (remaining cores should be a factor of $nWorkerCores)
nSchedulerCores=31
nClientCores=1
nWorkerCores=1	# should be a factor of 32 to ensure all available resources are utilized and no workers are spread across two compute nodes
WorkingDirName=testing

# set active directories
RUN_DIR=/path/to/$WorkingDirName
SCHEDULER_FILE=${RUN_DIR}/scheduler_file.json
SRC_DIR=/path/to/repo/

QUERY_LIST=/path/to/list_file.lst
TARGET_LIST=/path/to/pdb70_2022_03_19_structures.lst

date

PrimaryNode=$(echo $(scontrol show hostnames $SLURM_JOB_NODELIST) | awk '{print $1;}')

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ccs/home/davidsonrb/Apps/miniconda3_ANDES/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ccs/home/davidsonrb/Apps/miniconda3_ANDES/etc/profile.d/conda.sh" ]; then
        . "/ccs/home/davidsonrb/Apps/miniconda3_ANDES/etc/profile.d/conda.sh"
    else
        export PATH="/ccs/home/davidsonrb/Apps/miniconda3_ANDES/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# activate the conda environment
conda activate openmm_andes

if [ ! -d "$RUN_DIR" ]
then
    mkdir -p $RUN_DIR
fi
cd $RUN_DIR

# determine how many workers are able to be run on the Primary/Head Node
let x=$SLURM_CPUS_ON_NODE y=$nSchedulerCores z=$nClientCores nAvailCores=x-y-z	# subtracting the scheduler and client from the available cores; use this value for spinning up workers on the primary/head node
let x=$nAvailCores y=$nWorkerCores nWorkersPrimary=x/y	# returns an integer value, rounded down; if 0, then no workers will be spun up on primary node

# determine how many workers are able to run on the rest of the compute resources
let x=$SLURM_JOB_NUM_NODES y=1 nNodes=x-y		# remove the primary/head node from further math
let x=$nNodes y=$SLURM_CPUS_ON_NODE z=$nWorkerCores nWorkers=x*y/z	# calculating the total number of workers to be spun up by srun; excluding workers spun up on the primary/head node

let x=$nWorkersPrimary y=$nWorkers totalWorkers=x+y	# just for printing;

echo "################################################################################"
echo "Using python: " `which python3`
echo "PYTHONPATH: " $PYTHONPATH
echo "SRC_DIR: " $SRC_DIR
echo "Scheduler file:" $SCHEDULER_FILE
echo "Workflow running on" $SLURM_JOB_NUM_NODES "nodes"
echo "Primary/Head node is" $PrimaryNode
echo "Scheduler given" $nSchedulerCores "core(s) on the Primary/Head node"
echo "Client given" $nClientCores "core(s) on the Primary/Head node"
echo $totalWorkers "workers are spun up with" $nWorkerCores "core(s) per worker"
echo "################################################################################"

# gathering process ids for each step of the workflow.
dask_pids=""

##
## Start the dask scheduler 
##
echo "Spinning up the dask-scheduler"
srun -n 1 -N 1 -c $nSchedulerCores --cpu-bind=threads -w $PrimaryNode \
	dask-scheduler --interface ib0 --no-dashboard --no-show --scheduler-file $SCHEDULER_FILE  > dask_scheduler.out 2>&1 &
dask_pids="$dask_pids $!"

###
### Start the dask workers
###
# check to see if there's space on the primary/head node for workers; if there is, then submit some workers on that node
if [ $nWorkersPrimary -ne 0 ]; then
	echo "Spinning $nWorkersPrimary workers up on the head node, each getting $nWorkerCores cores"
	srun -n $nWorkersPrimary -N 1 -c $nWorkerCores --cpu-bind=threads -w $PrimaryNode \
		dask-worker --nthreads 1 --nworkers 1 --interface ib0 --no-dashboard --no-nanny --reconnect --scheduler-file ${SCHEDULER_FILE} > dask_worker.out 2>&1 &	
	dask_pids="$dask_pids $!"
fi

# check to see if there are other compute nodes in the allocation; if there are, then submit workers to those resources
if [ $nWorkers -ne 0 ]; then
	echo "Spinning $nWorkers workers up on non-head node(s), each getting $nWorkerCores cores"
	srun -n $nWorkers -N $nNodes -c $nWorkerCores --cpu-bind=threads \
		dask-worker --nthreads 1 --nworkers 1 --interface ib0 --no-dashboard --no-nanny --reconnect --scheduler-file ${SCHEDULER_FILE} > dask_worker.out 2>&1 &	
	dask_pids="$dask_pids $!"
fi

###
### Run the client task manager 
###
echo "Running the client script"
					  #--subdirectory-string USalign_ssi_lost_structures_pdb70_2022_03_19 \
					  #--subdirectory-string USalign_ssi_scope2.08_astral_rbd \
srun -n 1 -N 1 -c $nClientCores --cpu-bind=threads -w $PrimaryNode \
	python3 ${SRC_DIR}/dask_tskmgr.py --scheduler-file $SCHEDULER_FILE \
					  --query-pdb-list-file $QUERY_LIST \
					  --target-pdb-list-file $TARGET_LIST \
					  --script-path $SRC_DIR/runUSalign_sNS.sh \
					  --working-dir ${RUN_DIR} \
					  --subdirectory-string USalign_ssi_pdb70_2022_03_19 \
					  --timings-file timings.csv \
					  --tskmgr-log-file tskmgr.log \

# shutting down dask scheduler and worker commands
for pid in $dask_pids
do
        kill $pid
done

echo Run finished.

date

