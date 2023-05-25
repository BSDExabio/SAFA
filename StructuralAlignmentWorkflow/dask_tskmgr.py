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
                                     --script-path /path/to/dir/script.sh 
                                     
                                     OPTIONAL PARAMETERS:
                                     --ranking-metric [avgTMscore,maxTMscore,
                                                        TMscore1,TMscore2]
                                     --metric-cutoff 0.7
                                     
                                     
                                     # NOTE: add parameters to control saving files
                                     --working-dir /path/to/dir/ 
                                     --subdirectory-string DIRNAME
                                     
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

from USalignParser.usalign_parser import parse_usalign_file, threshold_comparison

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

def submit_pipeline(alignments, script, working_dir, scoring_metric, cutoff_threshold):
    """
    alignment task function for the dask workflow
    INPUT:
        :param alignments: list of shape Nx2, where N alignments occur between 
                  the two proteins in each N element
                  elem[0] is the path string pointing at the query protein
                  elem[1] is the path string associated with target protein,
                          being aligned to.
        :param script: path string pointing to the script that is executed in 
                       the subprocess call
        :param working_dir: string for path where temporary files 
                            will be written
        :param scoring_metric: string to be used as the scoring metric accepted
                               values = ''TMSCORE1', 'TMSCORE2', 'MAXTMSCORE',
                               'MINTMSCORE', 'AVGTMSCORE' (or any other combo of
                               capitalization)
        :param cutoff_threshold: float to be used as the cutoff for saving files
    RETURNS:
        :return: platform information i.e. "hostname"
        :return: worker identification string
        :return: start time, end time: epoch time values
        :return: number of keys in the results dict i.e. "nSuccesses"
        :return: the results_dict, dictionary with keys pointing to the target
                          and query protein, values being a list of quant
                          metrics; if an alignment passes the cutoff, this list
                          also includes the stdout from the alignment script
    """
    worker = get_worker()
    start_time = time.time()
    results_dict = {}

    # check to see if any scoring metric is considered
    if scoring_metric not in ['TMscore1','TMscore2','MaxTMscore','MinTMscore', 'AvgTMscore']:
        scoring_bool = False
    # check to see if a cutoff_threshold is defined
    elif not cutoff_threshold:
        scoring_bool = False
    else:
        scoring_bool = True
        cutoff_threshold = float(cutoff_threshold)

    for aln in alignments:
        # gather the paths to the structure files
        protein1 = aln[0]   # query
        protein2 = aln[1]   # target

        # create a temp string to save the trans/rot file under
        aln_uuid = working_dir + str(uuid4())
        # run the alignment script
        try:
            completed_process = subprocess.run(f'bash {script} {protein1} {protein2} {aln_uuid}',shell=True,capture_output=True,check=True)
        # if the alignment script errors out, skip the parsing steps
        except Exception as e:
            print(f'submit_pipeline, aln script, {protein1}, {protein2}', e, file=sys.stderr, flush=True)
            continue
        
        # parse the alignment output
        try:
            # creates a file-like-object text stream
            stdout_string = completed_process.stdout.decode() 
            stdout_file = io.StringIO(stdout_string)
            # put this text stream through the parser
            # hard coding the aln_algo as '' because we are not interested in 
            # the residue mapping at the present
            parsed_results = parse_usalign_file(stdout_file,'')[0]
            # save alignment results as a key-value in the results_dict
            # key: string with format of path/to/protein1 '|' path/to/protein2
            # value: list, [TMscore1, TMscore2, RMSD, SeqIDAli, Len1, Len2, 
            #               LenAligned]
            key = f'{protein1}|{protein2}'
            results_dict[key] = [float(parsed_results['TMscore1']),
                                 float(parsed_results['TMscore2']),
                                 float(parsed_results['RMSD']),
                                 float(parsed_results['SeqIDAli']),
                                 float(parsed_results['Len1']),
                                 float(parsed_results['Len2']),
                                 float(parsed_results['LenAligned'])]
        except Exception as e:
            print(f'submit_pipeline, parsing, {protein1}, {protein2}', e, file=sys.stderr, flush=True)
            continue
        
        # check to see if scoring should be considered; 
        # skip to next alignment if not True
        if not scoring_bool:
            continue
        # consider scoring_metric and cutoff_threshold to see if the 
        # stdout_string should be saved
        try:
            # check to see if the alignment's quantitative results pass the 
            # cutoff. if they do, then append the standard out string to the
            # results_dict[key] list. Otherwise, move on.
            cutoff_bool = threshold_comparison(parsed_results,ranking_metric,cutoff_threshold)
            if cutoff_bool:
                results_dict[key].append(stdout_string)
        except Exception as e:
            print(f'submit_pipeline, scoring, {protein1}, {protein2}', e, file=sys.stderr, flush=True)
            continue

    stop_time = time.time()
    return platform.node(), worker.id, start_time, stop_time, len(results_dict), results_dict


def target_identifier(target_str):
    """
    function to get an abbreviated string that will be used in the naming of 
    output directory paths
    """
    return Path(target_str).parent.name


def query_identifier(query_str):
    """
    function to get an abbreviated string that will be used in the naming of 
    output directory paths
    """
    return Path(query_str).stem


def post_analysis_pipeline(target_str, result_files_list, outputdir_str, subdir_str, sorted_by):
    """
    outputting combined results task function for the dask workflow
    INPUT:
        :param target_str: path string used to denote which target to collect
                           data for.
        :param results_files_list: list of path strings pointing to the results
                                   files that were created during the alignment
                                   task steps
        :param outputdir_str: path string for where results should be stored
        :param subdir_str: string for name used to create a subdirectory 
        :param sorted_by: string to be used as the scoring metric, accepted
                          values = ''TMSCORE1', 'TMSCORE2', 'MAXTMSCORE',
                          'MINTMSCORE', 'AVGTMSCORE' (or any other combo of 
                          capitalization)
    RETURNS:
        :return: start time, end time: epoch time values
        :return: string to denote which target was parsed during the task
    """
    start_time = time.time()
    # gather identifier for target
    target_id = target_identifier(target_str)
    # create subdirectory tree within which results associated with target_id 
    # will be written
    temp_path = Path(outputdir_str) / target_id / subdir_str
    temp_path.mkdir(mode=0o777,parents=True,exist_ok=True)
    # open and write alignment results to file
    with open(str(temp_path) + '/alignment_results.dat','w') as out_file, open(str(temp_path) + '/alignment_results.log','w') as log_file:
        out_file.write(f'Query,TMscore1,TMscore2,RMSD,SeqIDAli,Len1,Len2,LenAligned\n')
        # loop over and open the pickled results files
        for result_file in result_files_list:
            # load the pickle file
            with open(result_file,'rb') as infile:
                temp_results = pickle.load(infile)
            # temp_results is a dictionary with keys (query_path|target_path) 
            # and values, being a list of floats with organization of: 
            # [TMscore1, TMscore2, RMSD, SeqIDAli, Len1, Len2, LenAligned], an
            # eight element may be present that is the standard out from the 
            # alignment script
            for key, value in temp_results.items():
                # check to see if the target string in the dictionary key;
                # if not, skip to the next key, value pair
                if target_str not in key:
                    continue
                # check to see if the length of the value is 7
                elif len(value) == 7:
                    query_str = key.split('|')[0]
                    query_id = query_identifier(query_str)
                    out_file.write(f'{query_id},{value[0]},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]}\n')
                # check to see if the length of the value is 8; elem[7] is the
                # standard output from the alignment script; should include the
                # atomic pairs for alignments as well as the translation and 
                # rotation matrix to recreate the alignment
                elif len(value) == 8:
                    query_str = key.split('|')[0]
                    query_id = query_identifier(query_str)
                    out_file.write(f'{query_id},{value[0]},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]}\n')
                    output_path = temp_path / query_id / 'AlnResults.dat'
                    with open(str(output_path),'w') as result_file:
                        result_file.write(value[7])
                # if all else fail, something went wrong...
                else:
                    log_file.write( f"{key}, {value}, len(value) does not match expected length of results list. Something wrong with the USalign results?\n")
                    continue
        
    # read in the results file and parse
    df = pandas.read_csv(str(temp_path) + '/alignment_results.dat')
    # create a maxTMscore column if expecting to sort by this value
    if sorted_by.upper() == 'MAXTMSCORE':
        df['maxTMscore'] = np.max(df[['TMscore1','TMscore2']],axis=1)
    # create a minTMscore column if expecting to sort by this value
    elif sorted_by.upper() == 'MINTMSCORE':
        df['maxTMscore'] = np.min(df[['TMscore1','TMscore2']],axis=1)
    # create a avgTMscore column if expecting to sort by this value
    elif sorted_by.upper() == 'AVGTMSCORE':
        df['avgTMscore'] = np.mean(df[['TMscore1','TMscore2']],axis=1)
    # check that sorted_by is one of the column names of the panda dataframe
    elif: #sorted_by not in list(df.columns):
        print(f'{sort_by} not in pandas dataframe columns {list(df.columns)}. No ranked_alignment_results.dat file is written.', file=sys.stderr, flush=True)
        stop_time = time.time()
        return start_time, stop_time, target_str

    # sort the dataframe
    sorted_df = df.sort_values(sorted_by,ascending=False)
    # save to file
    sorted_df.to_csv(str(temp_path) + '/ranked_alignment_results.dat',index=False)

    stop_time = time.time()
    return start_time, stop_time, target_str


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Molecular dynamics simulation task manager')
    # essential parameters
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    parser.add_argument('--query-pdb-list-file', '-inp', required=True, help='list file that contains the paths and sequence lengths of protein models that are aligned to the target structures; the structural library')
    parser.add_argument('--target-pdb-list-file', '-lib', required=True, help='list file that contains the paths and sequence lengths of protein models associated with targets being aligned to')
    parser.add_argument('--script-path', '-sp', required=True, help='path that points to the script for the subprocess call')
    parser.add_argument('--subdirectory-string', '-sub', required=True, help='string to be used for consistent naming of subdirectories within which results files will be saved for each target structure.')
    # optional parameters
    parser.add_argument('--working-dir', '-wd', required=False, default='./', help="path that points to the working directory for the output files; default = './'")
    parser.add_argument('--scoring-metric', '-metric', required=False, default=None, help='string used to denote which alignment column to rank structures by; default is empty')
    parser.add_argument('--cutoff-threshold', '-thresh', required=False, default=None, help='float value to be used to determine which alignments have high-verbosity results written to storage.')
    args = parser.parse_args()

    # check if working_dir is already made
    if not Path(args.working_dir).is_dir():
        Path(args.working_dir).mkdir(parents=True)

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Query models are listed in {args.query_pdb_list_file}')
    main_logger.info(f'Target models are listed in {args.target_pdb_list_file}')
    main_logger.info(f'Path to subprocess script: {args.script_path}')
    main_logger.info(f'Working directory: {args.working_dir}')
    main_logger.info(f'Alignment results will be saved within subdirectories named {args.subdirectory_string}')
    if args.scoring_metric:
        main_logger.info(f'All alignment results will be ranked by the {args.scoring_metric} metric')
    if args.cutoff_threshold:
        main_logger.info(f'The metric threshold value of {args.cutoff_threshold} will be used to save high-quality alignment results to storage.')
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

    # parsing the query_pdb_list_file
    # query_list represents the list of structures that will be aligned to the 
    # target structures
    main_logger.info(f'Reading the query structure file.')
    with open(args.query_pdb_list_file,'r') as structures_file:
        query_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    
    # sorted largest to smallest of the structure list; system size is a basic but good estimate of computational cost
    sorted_query_list = sorted(query_list, key = lambda x: int(x[1]))[::-1]
    del query_list
    sorted_query_list = [elem[0] for elem in sorted_query_list]
    main_logger.info(f'Preparing to run {len(sorted_query_list):,} query structures, sorted by longest to shortest sequences.')

    # parse the target_pdb_list_file
    # target_list represents the list of structures that query structures will 
    # align to
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
    main_logger.info(f'A total of {nAlignments:,} alignments will be performed.')

    ############
    # if total number of alignment tasks is relatively small, 
    # then create the full iterable list and map directly to client
    ############
    if nAlignments < 10**4:
        main_logger.info(f'The total number of alignments is relatively small. Running all alignments as individual tasks.')
        
        # do the thing.
        alignments_list = list(itertools.product(sorted_query_list,sorted_target_list))
        aln_futures = client.map(submit_pipeline, alignments_list, script = args.script_path, working_dir = args.working_dir, scoring_metric = args.scoring_metric, cutoff_threshold = args.cutoff_threshold, pure=False)

        # gather alignment results.
        results_dict = {}
        aln_completed = as_completed(aln_futures)
        for finished_task in aln_completed:
            hostname, worker_id, start_time, stop_time, nSuccesses, results = finished_task.result()
            query, target = list(results.keys())[0].split('|')
            main_logger.info(f'Alignment between {query} and {target} has finished.')
            append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, 1)
            results_dict.update(results)
        
        out = f'{args.working_dir}/results_dictionary.pkl'
        main_logger.info(f'RESULTS WRITTEN TO {out}. {time.time()}')
        with open(out,'wb') as out_file:
            pickle.dump(results_dict,out_file)

        parse_futures = client.map(post_analysis_pipeline, sorted_target_list, result_files_list = [out], outputdir_str = args.working_dir, subdir_str = args.subdirectory_string, sorted_by = args.scoring_metric, pure=False)
        parse_completed = as_completed(parse_futures)
        for finished_task in parse_completed:
            start_time, stop_time, target = finished_task.result()
            main_logger.info(f'Parsed alignment hits for {target}, taking {stop_time - start_time} seconds.')

    ############
    # if total number of tasks is not relatively small, 
    # then partition the query list and map query-target partitions to the client
    ############
    else:
        main_logger.info(f'The total number of alignments ({nAlignments:,}) is large. Will run the alignments in batches of query structures.')
        # get number of workers associated with client
        NUM_WORKERS = len(client.scheduler_info()['workers'].keys())
        # break the query_list into sublists
        query_sublists = [sorted_query_list[i::NUM_WORKERS*5] for i in range(NUM_WORKERS*5)]
        # approx how many structures in each sublist
        nQueries_per_sublist = int(np.mean([len(sublist) for sublist in target_sublists]))
        # calc the number of tasks, where tasks are now multiple alignments 
        nTasks = len(sorted_target_list)*len(query_sublists)
        main_logger.info(f'The total number of tasks is {nTasks:,} where each task is the alignment of multiple query structures to target structures. Average lengths of sublists is {nQueries_per_sublist}).')
        # loop over target structures; submit tasks associated with that target
        for target in sorted_target_list:
            target_start = time.time()

            # do the thing.
            alignment_lists = [list(zip(query_sublist,[target]*len(query_sublists))) for query_sublist in query_sublists]
            aln_futures = client.map(submit_pipeline, alignment_list, script = args.script_path, working_dir = args.working_dir, scoring_metric = args.scoring_metric, cutoff_threshold = args.cutoff_threshold, pure=False)
            
            # gather target's alignment results.
            results_dict = {}
            aln_completed = as_completed(aln_futures)
            for finished_task in aln_completed:
                hostname, worker_id, start_time, stop_time, nSuccesses, results = finished_task.result()
                append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, nSuccesses)
                results_dict.update(results)
        
            out = f'{args.working_dir}/results_dictionary.pkl'
            with open(out,'wb') as out_file:
                pickle.dump(results_dict,out_file)
    
            post_analysis_pipeline(target, result_files_list = [out], outputdir_str = args.working_dir, subdir_str = args.subdirectory_string, sorted_by = args.scoring_metric, pure=False)
            main_logger.info(f'Finished running, collecting, and outputting alignment results for {target}. This took {time.time() - target_start} seconds.')
    
    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

