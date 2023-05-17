#!/usr/bin/env python3
""" Task manager for running the dask pipeline for running structure alignment 
    of a set of target protein structures against a library of pdb structures.
    Save the alignments if so desired.
    USAGE: 
        python3 dask_taskmgr.py [-h] --scheduler-file SCHEDULER_FILE 
                                     --target-pdb-list-file INPUT_FILE 
                                     --query-pdb-list-file LIB_DIR_PATH
                                     --timings-file TIMINGS_FILE.csv 
                                     --tskmgr-log-file TSKMGR.log
                                     --ranking-metric [avgTMscore,maxTMscore,
                                                        TMscore1,TMscore2]
                                     --metric-cutoff 0.7
                                     --script-path /path/to/dir/script.sh 
                                     
                                     --num-workers nWorkers
                                     --working-dir /path/to/dir/ 
                                     --subdirectory-string DIRNAME
                                     # NOTE: add parameters to control saving files
    INPUT: 
        -h, --help      show this help message and exit
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
        --target-pdb-list-file INPUT_FILE, -inp INPUT_FILE
                        list file that contains the paths and sequence lengths 
                        of protein models
        --query-pdb-list-file LIB_DIR_PATH, -lib LIB_DIR_PATH
                        path to a directory within which a library of protein 
                        models are held
        --timings-file TIMINGS_FILE.csv, -ts TIMINGS_FILE.csv 
                        CSV file for task processing timings
        --tskmgr-log-file TSKMGR.log, -log TSKMGR.log 
                        path for a file within which logging output for the 
                        workflow will be written
        --ranking-metric COLUMN_NAME, -rank COLUMN_NAME
                        string to set which column of the pandas dataframe to 
                        rank the alignment hits with. options: 'avgTMscore', 
                                                               'maxTMscore', 
                                                               'TMscore1', 
                                                               'TMscore2'
        --metric-cutoff FLOAT_VALUE, -cutoff FLOAT_VALUE
                        float value set to determine the ranking metric value 
                        above which, alignment results will be saved.
        --script-path /path/to/dir/script.sh, -sp /path/to/dir/script.sh
                        full path to the script to be run within the subprocess
       


        --working-dir /path/to/dir/, -wd /path/to/dir/
                        full path to the directory within which files will be 
                        written
        --subdirectory-string DIRNAME, -sub DIRNAME
                        string to be used for consistent naming of subdirectories 
                        within which results files will be saved for each query
                        structure.
"""

import sys
import time
import pickle
import argparse
import platform
import logging
from uuid import uuid4
import itertools
import numpy as np
from pathlib import Path
import pandas
import csv

import subprocess
from subprocess import CalledProcessError

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

#######################################
### LOGGING FUNCTIONS
#######################################

def append_timings(csv_writer, file_object, hostname, worker_id, start_time, 
                   stop_time, nSuccesses):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param nSuccesses: symbolic result of the subprocess call
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'nSuccesses' : nSuccesses})
    file_object.flush()


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


#######################################
### DASK RELATED FUNCTIONS
#######################################

def submit_pipeline(proteins, script, ranking_metric, cutoff_threshold, working_dir='./'):
    """
    alignment task function for the dask workflow
    INPUT:
        :param proteins: list with length of 2
                  proteins[0] is the path string pointing at the target protein
                  proteins[1] is a path string or a list of path strings 
                              associated with query proteins
        :param script: path string pointing to the script that is executed in 
                       the subprocess call
        :param ranking_metric: string to be used with the cutoff
        :param cutoff_threshold: float to be used as the cutoff for saving files
        :param working_dir: string for path where temporary files will be written
    RETURNS:
        :return: platform information i.e. "hostname"
        :return: worker identification string
        :return: start time, end time: epoch time values
        :return: number of keys in the results dict i.e. "nSuccesses"
        :return: the results_dict, dictionary with keys pointing to the target
                          and queryproteins, values being a list of quant
                          metrics if an alnment passes the cutoff, this list
                          also includes the stdout from the script
    """
    worker = get_worker()
    start_time = time.time()
    results_dict = {}
    target_protein = proteins[0]
    query_proteins = proteins[1]
    if type(query_proteins) == str:
        query_proteins = [query_proteins]

    for query in query_proteins:
        # create a temp string to save the trans/rot file under
        aln_uuid = working_dir + str(uuid4())
        try:
            completed_process = subprocess.run(f'bash {script} {query} {target_protein} {aln_uuid}',shell=True,capture_output=True,check=True)
            # creates a file-like-object text stream
            stdout_string = completed_process.stdout.decode() 
            stdout_file = io.StringIO(stdout_string)
            # put this text stream through the parser
            # hard coding the aln_algo as '' because we are not interested in 
            # the residue mapping at the present
            parsed_results = usalign_parser.parse_usalign_file(stdout_file,'')[0]
            # save alignment results as a key-value in the results_dict
            # key: string with format of path/to/target/protein '_' path/to/query/protein
            # value: list, [TMscore1, TMscore2, RMSD, SeqID1, SeqID2, AlnSeqID, Len1, Len2, alnLen]
            results_dict[f'{target_protein}_{query}'] = [float(parsed_results['tmscore1']),
                                                         float(parsed_results['tmscore2']),
                                                         float(parsed_results['rmsd']),
                                                         float(parsed_results['alnSeqID']),
                                                         float(parsed_results['Len1']),
                                                         float(parsed_results['Len2']),
                                                         float(parsed_results['alnLen'])]
            # check to see if the alignment's quantitative results pass the 
            # cutoff. if they do, then append the standard out string to the
            # results_dict[key] list. Otherwise, move on.
            cutoff_bool = usalign_parser.threshold_comparison(parsed_results,ranking_metric,cutoff_threshold)
            if cutoff_bool:
                results_dict[f'{target_protein}_{query}'].append(stdout_string)

        except CalledProcessError as e:
            print(query_protein, target, e, file=sys.stderr, flush=True)
        
    stop_time = time.time()
    return platform.node(), worker.id, start_time, stop_time, len(results_dict), results_dict


def post_analysis_pipeline(query_str, result_files_list, outputdir_str= './', subdir_str='TMalign', nRanked=100, sorted_by='maxTMscore'):
    """
    """
    start_time = time.time()
    protein_id = Path(query_str).parent.name
    temp_path = Path(outputdir_str) / protein_id / subdir_str
    temp_path.mkdir(mode=0o777,parents=True,exist_ok=True)
    # open and write alignment results to file
    with open(str(temp_path) + '/alignment_results.dat','w') as out_file, open(str(temp_path) + '/alignment_results.log','w') as log_file:
        out_file.write(f'Target Path,TMscore1,TMscore2,RMSD,SeqIDali,Len1,Len2,LenAligned\n')
        for result_file in result_files_list:
            with open(result_file,'rb') as infile:
                temp_results = pickle.load(infile)
            for key, value in temp_results.items():
                if query_str not in key:
                    continue
                elif len(value) != 9:
                    log_file.write( f"{key}, {value}, len(value) does not match expected length. something wrong with the TMalign results?\n")
                    continue
                else:
                    target = '/' + key.split('_/')[1]
                    out_file.write(f'{target},{value[0]},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]},{value[7]},{value[8]}\n')
        
    # read in the results file and parse
    df = pandas.read_csv(str(temp_path) + '/alignment_results.dat')
    # create a maxTMscore column if expecting to sort by this value
    if sorted_by == 'maxTMscore':
        df['maxTMscore'] = np.max(df[['TMscore1','TMscore2']],axis=1)
    elif sorted_by == 'avgTMscore':
        df['avgTMscore'] = np.mean(df[['TMscore1','TMscore2']],axis=1)
    # check that sorted_by is one of the column names of the panda dataframe
    elif sorted_by not in list(df.columns):
        print(f'{sort_by} not in pandas dataframe columns {list(df.columns)}. No ranked_alignment_results.dat file is written.', file=sys.stderr, flush=True)
        stop_time = time.time()
        return start_time, stop_time, query_str

    # sort the dataframe
    sorted_df = df.sort_values(sorted_by,ascending=False)
    # write the sorted dataframe to file
    try:
        sorted_df.head(nRanked).to_csv(str(temp_path) + '/ranked_alignment_results.dat',index=False)
    # if nRanked > the total number of alignments performed, the above code will raise an exception
    # instead, just write the full set of ranked alignment results. 
    except:
        sorted_df.to_csv(str(temp_path) + '/ranked_alignment_results.dat',index=False)

    stop_time = time.time()
    return start_time, stop_time, query_str


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Molecular dynamics simulation task manager')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--query-pdb-list-file', '-inp', required=True, help='list file that contains the paths and sequence lengths of protein models')
    parser.add_argument('--target-pdb-list-file', '-lib', required=True, help='list file that contains the paths and sequence lengths of protein models assocaited with the target library')
    parser.add_argument('--script-path', '-sp', required=True, help='path that points to the script for the subprocess call')
    parser.add_argument('--working-dir', '-wd', required=True, help='path that points to the working directory for the output files')
    parser.add_argument('--subdirectory-string', '-sub', required=True, help='string to be used for consistent naming of subdirectories within which results files will be saved for each query structure.')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    parser.add_argument('--nRanked-structures', '-nrank', required=True, help='integer value set for number of top scoring alignment hits to ranked file')
    parser.add_argument('--sorted-by', '-sort', required=True, help='string used to denote which column to rank structures by')
    parser.add_argument('--num-workers', '-nw', required=True, help='integer used to set how many workers are being claimed by the client')
    args = parser.parse_args()

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Working directory: {args.working_dir}')
    main_logger.info(f'Path to subprocess script: {args.script_path}')
    main_logger.info(f'Query models are listed in {args.query_pdb_list_file}')
    main_logger.info(f'Target models are listed in {args.target_pdb_list_file}')
    main_logger.info(f'Alignment results will be saved within subdirectories named {args.subdirectory_string}')
    main_logger.info(f'All alignment results will be ranked by the {args.sorted_by} metric')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')
    
    # set up timing log file.
    main_logger.info(f'Opening the timing file.')
    timings_file = open(args.timings_file, 'w')
    timings_csv  = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','nSuccesses'])
    timings_csv.writeheader()

    # parse the query_pdb_list_file
    main_logger.info(f'Reading the query structure file.')
    with open(args.query_pdb_list_file,'r') as structures_file:
        query_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    
    # sorted largest to smallest of the structure list; system size is a basic but good estimate of computational cost
    sorted_query_list = sorted(query_list, key = lambda x: int(x[1]))[::-1]
    del query_list
    sorted_query_list = [elem[0] for elem in sorted_query_list]
    main_logger.info(f'Preparing to run {len(sorted_query_list):,} query structures, sorted by longest to shortest sequences.')

    # parse the target_pdb_list_file
    main_logger.info(f'Reading the target structure file.')
    with open(args.target_pdb_list_file,'r') as structures_file:
        target_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    
    # sorted largest to smallest of the structure list; system size is a basic but good estimate of computational cost
    sorted_target_list = sorted(target_list, key = lambda x: int(x[1]))[::-1]
    del target_list
    sorted_target_list = [elem[0] for elem in sorted_target_list]
    nTargets = len(sorted_target_list)
    main_logger.info(f'Preparing to run alignments against {nTargets:,} target structures.')
    
    # determining total number of alignments to be perforemd
    nAlignments = len(sorted_query_list)*len(sorted_target_list)
    main_logger.info(f'A total of {nAlignments} alignments will be performed.')
   




    ############
    # if total number of alignment tasks is relatively small, 
    # then create the full iterable list and map directly to client
    ############
    if nAlignments < 10**4 or nTargets < NUM_WORKERS:
        main_logger.info(f'The total number of alignments ({nAlignments:,}) is relatively small. Running all alignments as individual tasks.')
        
        # do the thing.
        aln_futures = client.map(submit_pipeline, list(itertools.product(sorted_query_list,sorted_target_list)), script = args.script_path, pure=False) #, batch_size=10**4

        # gather alignment results.
        results_dict = {}
        aln_completed = as_completed(aln_futures)
        for finished_task in aln_completed:
            hostname, worker_id, start_time, stop_time, nSuccesses, results = finished_task.result()
            query_target = list(results.keys())[0]
            main_logger.info(f'{query_target} alignment has finished.')
            append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, query_target)
            results_dict.update(results)
        
        out = f'{args.working_dir}/results_dictionary.pkl'
        main_logger.info(f'RESULTS WRITTEN TO {out}. {time.time()}')
        with open(out,'wb') as out_file:
            pickle.dump(results_dict,out_file)

        parse_futures = client.map(post_analysis_pipeline, sorted_query_list, result_files_list = [out], outputdir_str = args.working_dir, subdir_str = args.subdirectory_string, nRanked = int(args.nRanked_structures), sorted_by = args.sorted_by, pure=False)
        parse_completed = as_completed(parse_futures)
        for finished_task in parse_completed:
            start_time, stop_time, query = finished_task.result()
            main_logger.info(f'Parsed alignment hits for {query}, taking {stop_time - start_time} seconds.')

    ############
    # if total number of tasks is not relatively small, 
    # then partition the target list and map query-target partitions to the client
    ############
    else:
        main_logger.info(f'The total number of alignments ({nAlignments:,}) is large. Will run the alignments in batches of targets.')
        target_sublists = [sorted_target_list[i::NUM_WORKERS*3] for i in range(NUM_WORKERS*3)]
        nTargets_per_sublist = int(np.mean([len(sublist) for sublist in target_sublists]))
        #nTargets_per_sublist = nTargets//NUM_WORKERS + (nTargets % NUM_WORKERS > 0)
        #target_sublists = [target_list[i:i+nTargets_per_sublist] for i in range(0,len(target_list),nTargets_per_sublist)] # not desired if target_list is sorted largest to smallest
        nTasks = len(sorted_query_list)*len(target_sublists)
        main_logger.info(f'The total number of tasks is {nTasks:,} where each task is the alignment of a query structure to a set of target structures (average lengths of {nTargets_per_sublist}).')
        for query in sorted_query_list:
            query_start = time.time()
            # do the thing.
            aln_futures = client.map(submit_pipeline, list(itertools.product([query],target_sublists)), script = args.script_path, pure=False) #, batch_size=10**4
            
            # gather query's alignment results.
            results_dict = {}
            aln_completed = as_completed(aln_futures)
            for finished_task in aln_completed:
                hostname, worker_id, start_time, stop_time, nSuccesses, results = finished_task.result()
                append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, nSuccesses)
                results_dict.update(results)
        
            out = f'{args.working_dir}/results.pkl'
            with open(out,'wb') as out_file:
                pickle.dump(results_dict,out_file)
    
            post_analysis_pipeline(query, result_files_list = [out], outputdir_str = args.working_dir, subdir_str = args.subdirectory_string, nRanked = int(args.nRanked_structures), sorted_by = args.sorted_by)
            main_logger.info(f'Finished running, collecting, and outputting alignment results for {query}. This took {time.time() - query_start} seconds.')
    
    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)
    #client.shutdown()
