
date

WorkingDirName=example

# set active directories
RUN_DIR=./$WorkingDirName
SRC_DIR=.

export PYTHONPATH=$PYTHONPATH:$SRC_DIR
blastp_command="/home/russ/Apps/ncbi-blast-2.13.0+/bin/blastp"

if [ ! -d "$RUN_DIR" ]
then
    mkdir -p $RUN_DIR
fi
#cd $RUN_DIR

echo "################################################################################"
echo "Using python: " `which python3`
echo "PYTHONPATH: " $PYTHONPATH
echo "SRC_DIR: " $SRC_DIR
echo "Rerun using: python3 ${SRC_DIR}/processing_pdb_library.py --input-list-file $RUN_DIR/example.lst --output-directory $RUN_DIR/example_output --timings-file timings.csv --tskmgr-log-file tskmgr.log --blastp-path $blastp_command"
echo "################################################################################"

###
### Run the client task manager 
###
echo "Running the client script"
python3 ${SRC_DIR}/processing_pdb_library.py --input-list-file $RUN_DIR/example.lst \
					     --output-directory $RUN_DIR/example_output \
					     --timings-file timings.csv \
					     --tskmgr-log-file tskmgr.log \
					     --blastp-path $blastp_command

echo Run finished.

date

