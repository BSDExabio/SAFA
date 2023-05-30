
date

# set active directories
RUN_DIR=./
SRC_DIR=..
WorkingDirName=large_test

QUERY_LIST=queries.lst
TARGET_LIST=large_targets.lst

if [ ! -d "$RUN_DIR" ]
then
    mkdir -p $RUN_DIR
fi

echo "################################################################################"
echo "Using python: " `which python3`
echo "PYTHONPATH: " $PYTHONPATH
echo "SRC_DIR: " $SRC_DIR
echo "Rerun using: python3 ${SRC_DIR}/dask_tskmgr.py "
echo "################################################################################"

###
### Run the client task manager 
###
echo "Running the client script"
python3 ${SRC_DIR}/dask_tskmgr.py --timings-file $WorkingDirName/timings.csv \
				  --tskmgr-log-file $WorkingDirName/tskmgr.log \
				  --query-pdb-list-file $QUERY_LIST \
				  --target-pdb-list-file $TARGET_LIST \
				  --script-path $SRC_DIR/runUSalign_sNS.sh \
				  --subdirectory-string sample_structure_library \
				  --working-dir ${WorkingDirName}/ \
				  --scoring-metric avgTMscore \
				  --cutoff-threshold 0.5

echo Run finished.

date

