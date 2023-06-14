
import sys
import time
import datetime
import pickle
import argparse
import platform
import logging
import warnings
from uuid import uuid4
from pathlib import Path
import pandas
import csv
from io import StringIO
import copy
import os

import MDAnalysis
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

import rcsb_query
import uniprot_query

#######################################
### LOGGING FUNCTIONS
#######################################

def append_timings(csv_writer, file_object, hostname, worker_id, start_time, 
                   stop_time, task_type, return_code):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param file_object: file object associated with the CSV
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param task_type: integer used to denote which step of the workflow has 
                      been performed
    :param return_code: result of the task; 1=successful, 0=failure
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'task_type'  : task_type,
                         'return_code': return_code})
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

def get_num_workers(client):
    """ Get the number of active workers
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())


#######################################
### SUB-FUNCTIONS
#######################################

def get_pdbid_chainid(pdb_file):
    """Gather the PDBID_CHAINID string from a pdb_file
    NOTE: this function might need to change based on where the PDBID_CHAINID
    information is stored within the pdb_file.
    
    :param pdb_file: path object that points to a pdb file, which has an 
                     associated PDBID_CHAINID.
    :return: string of format AAAA_B where AAAA is the four character RCSB PDB
             ID and B is the chain ID of interest
    
    For the PDB70 structure set, pdb files have the PBDID_CHAINID as the stem.
    For the SCOPe structure set, PDBID_CHAINID are stored in the REMARK lines.
    """
    # FOR PARSING PDB70 structures
    return pdb_file.stem
    # FOR PARSING SCOPe structures
    #return ....

#######################################
### TASK FUNCTIONS
#######################################
# Since tasks 0 and 1 are gonna be strung together and parsed in the same 
# `as_completed' for-loop, we need their function's output to be 
# similarly formatted, even if these similarities are superficial.
# Output format: 
#       task_num: a number to denote which task in the workflow has been 
#                 completed.
#       results: a dictionary. The inner dictionary can be filled with any 
#                key:value pairs; the key(s) for the outer dictionary should
#                be the relevant ID(s) that was parsed.
#       hostname: a string that describes the platform/compute node the task
#                 was performed on.
#       worker_id: a string that describes the dask worker that performed
#                  the task
#       start,stop: epoch times for when the task began and finished.
#       return_code: integer value, values dependent by function.

# step0
def _gather_pdb_file(PDBID_chainIDs, output_directory = './'):
    """ download the structure files and gather uniprotkb accession id
    :param PDBID_chainIDs: list, formatted as [string, list of strings].
                           elem[0] is a PDB accession ID. 
                           elem[1] is a list of strings.
    :param output_directory: string, local or global path string to a directory
                             within which structural files will be saved.
    :return: step_number, subdictionary, host_information, worker_id,
             start_time, stop_time, n_successes
             
             step_number: int, just denoting which step in the process has been
                          completed; 0 for the _parse_pdb_file function. 
             subdictionary: dictionary of a dictionary; format:
                            {'PDBID_CHAINID': {'structure_path': str,
                                               'struct_first_resid': int,
                                               'struct_last_resid': int,
                                               'struct_seq': str,
                                               'UniProtID': str}
                             ...}
             host_information: string, denotes what hardware resources were 
                               used for this task.
             worker_id: string, denotes which dask worker performed this task.
             start_time,stop_time: floats, epoch times when the task began and 
                                   finished, respectively.
             n_successes: int, nChainIDs if all chains are successfully 
                          gathered and parsed, < nChainIDs otherwise
    """
    start_time = time.time()
    n_successes = 0
    worker = get_worker()
    
    # parse input list
    PDBID    = PDBID_chainIDs[0]    # string of RCSB accession ID
    chainIDs = PDBID_chainIDs[1]    # list of strings, chain IDs

    # check to see if the output_directory already exists, if not, then break 
    # out of this task; this will force users to create the directory before
    # submitting this type of task
    if not Path(output_directory).is_dir():
        print(f'{output_directory} does not exist, cannot write structure files. Exception: {e}', file=sys.stderr, flush=True)
        return 0, {PDBID:'failed to be gathered'}, platform.node(), worker.id, start_time, time.time(), n_successes

    # gather the originating PDB structure, all chains included
    mda_universe = rcsb_query.gather_structure(PDBID)
    # if mda_universe is None, then no structure was gathered, return an empty 
    # dictionary result
    if not mda_universe:
        return 0, {PDBID:'failed to be gathered'}, platform.node(), worker.id, start_time, time.time(), n_successes

    # else, start to parse the chainIDs and gather non-structural data
    pdb_dict = {}
    file_names = []

    # setting chainID attribute to the MDA universe object since that data
    # structure does not include that metadata by default. Instead, MDA uses
    # the `segid' classifier/atom selection string.
    # We're cheating a bit here... setting all atoms in the universe to have
    # `chainID' as the topology attribute chainID; setting to '0' temporarily
    mda_universe.add_TopologyAttr('chainID',len(mda_universe.atoms)*['0'])
    
    # loop over the chainIDs of interest and gather the structures and metadata
    for i, chainID in enumerate(chainIDs):
        # setting the structure_id that will be used for the remainder of the 
        # dataset; format will be RCSBAccessionID_ChainID string
        structure_id = f'{PDBID}_{chainID}'
        
        # gather the structure's associated uniprotkb accession id. If no id
        # is gathered or the function fails, the structure's uniprotkb id will
        # be None. 
        try:
            uniprotid = rcsb_query.query_uniprot_str(structure_id)
        except Exception as e:
            print(f'failed to pull the UniProt accession id associated with {structure_id}. Exception: {e}', file=sys.stderr, flush=True)
            uniprotid = None

        # change the chainID string to the appropriate one
        mda_universe.atoms.chainIDs = len(mda_universe.atoms)*[chainID]

        # Making the atom selection that will be written to file. Need to
        # check that only one model is written to file. NMR structures have 
        # numerous models that are not considered "frames" in a trajectory when 
        # loaded into MDAnalysis. Only grabbing the first model out of any set 
        # of structural models. 
        if len(mda_universe.models) > 1:
            # only append the remark about grabbing the first model's structure 
            # once
            if i == 0:
                mda_universe.trajectory.remarks.append(f'NOTE: only grabbing the atoms of the 1st model associated with this PDBID')
            sel = mda_universe.models[0].select_atoms(f'protein and segid {chainID} and not resname ACE ALAD ASF CME CYX DAB NME QLN ORN')
        # otherwise, we can just make the desired atom selection
        else:
            sel = mda_universe.select_atoms(f'protein and segid {chainID} and not resname ACE ALAD ASF CME CYX DAB NME QLN ORN')
        
        # add a remark about the chain to be written to file, will always be the
        # last remark
        if i == 0:
            mda_universe.trajectory.remarks.append(f'NOTE: only grabbing atoms associated with chain/segment {chainID}')
            mda_universe.trajectory.remarks.append(f'UniProtKB Accession ID: {uniprotid}')
        # replace the already present remark
        else:
            mda_universe.trajectory.remarks[-2] = f'NOTE: only grabbing atoms associated with chain/segment {chainID}'
            mda_universe.trajectory.remarks[-1] = f'UniProtKB Accession ID: {uniprotid}'
        
        # write structure to file
        output_file = f'{output_directory}{PDBID}_{chainID}.pdb'
        with warnings.catch_warnings():
            # ignore some annoying warnings from sel.write line due to missing information (chainIDs, elements, and record_types). 
            warnings.simplefilter('ignore',UserWarning)

            sel.write(f'{output_file}',bonds=None)
        
        # fill out the pdb_dict object with the structure_id's metadata
        pdb_dict[structure_id] = {}
        pdb_dict[structure_id]['structure_path'] = output_file
        pdb_dict[structure_id]['struct_first_resid'] = sel.residues[0].resid
        pdb_dict[structure_id]['struct_last_resid']  = sel.residues[-1].resid
        pdb_dict[structure_id]['struct_seq']  = sel.residues.sequence(format='string')
        pdb_dict[structure_id]['UniProtID']  = uniprotid

        n_successes += 1

    return 0, pdb_dict, platform.node(), worker.id, start_time, time.time(), n_successes


# step 1
def _query_uniprot_flat_file(uniprotid):
    """Gather UniProt flat file metadata associated with a uniprotid
    :param uniprotid: string, numerous formats; expected to be associated with 
                      an entry in the UniProtKB database
    :return: step_number, subdictionary, host_information, worker_id, 
                start_time, stop_time, return_code
             
             step_number: int, just denoting which step in the process has been
                          completed; 2 for the _query_uniprot_flat_file func. 
             subdictionary: dictionary of a dictionary; format:
                            {'UniProtID': meta_dict}
                                meta_dict is of a very complex form. See 
                                uniprot_query.request_uniprot_metadata for 
                                further details
             host_information: string, denotes what hardware resources were 
                               used for this task.
             worker_id: string, denotes which dask worker performed this task.
             start_time,stop_time: floats, epoch times when the task began and 
                                   finished, respectively.
             return_code: int, 0 if successfully completed, 1 or -9999 otherwise
    """
    start_time = time.time()
    worker = get_worker()
    meta_dict = {}
    return_code = -9999
    try:
        meta_dict = uniprot_query.request_uniprot_metadata(uniprotid)
        if meta_dict['query_status'] == 200:
            return_code = 0
        else: 
            return_code = 1
    except Exception as e:
        print(f'failed to pull {uniprotid} flat file. Exception: {e}', file=sys.stderr, flush=True)

    return 1, {uniprotid: meta_dict}, platform.node(), worker.id, start_time, time.time(), return_code


# step 2
def _blast_aln(structure_seq,uniprotid_seq,blastp_path,ids):
    """Run and parse a BlastP sequence alignment between the structure and flat 
       file sequences

    :param structure_seq: str, sequence from the PDBID_CHAINID structure 
    :param uniprotid_seq: str, sequence from the UniProtID flat file
    :param blastp_path: string, global path to the BlastP executable
    :param ids: list of strings, format: [PDBID_CHAINID,UniProtID]
    
    :return: step_number, subdictionary, ids, host_information, worker_id, 
                start_time, stop_time, return_code
             
             step_number: int, just denoting which step in the process has been
                          completed; 2 for the _blast_aln func. 
             subdictionary: dictionary; format:
                            {'query_start': int,
                             'query_end': int,
                             'sbjct_start': int,
                             'sbjct_end': int,
                             'alignment': list of tuples of len 3, format: 
                                          (query_resname, match character, 
                                           sbjct resname),
                             'e-value': float,
                             'score': float,
                             'bits': float,
                             'n_gaps': int}
             ids: list of strings, format: [PDBID_CHAINID,UniProtID]
             host_information: string, denotes what hardware resources were 
                               used for this task.
             worker_id: string, denotes which dask worker performed this task.
             start_time,stop_time: floats, epoch times when the task began and 
                                   finished, respectively.
             return_code: int, 0 if successfully completed, -9999 otherwise
    """
    start_time = time.time()
    worker = get_worker()
    seq_aln_dict = {}
    return_code = -9999
    
    # unfortunately need to write fasta files for blastp to run
    # query will always be the PDBID_CHAINID structure's sequence
    # subject will always be the uniprot flat file's sequence
    # if files cannot be written, then this task will return
    try:
        aln_uuid = str(uuid4())
        with open(f'{aln_uuid}_query.fasta','w') as fsta:
            fsta.write(f">query\n{structure_seq}")
        with open(f'{aln_uuid}_sbjct.fasta','w') as fsta:
            fsta.write(f"{uniprotid_seq}")
    except Exception as e:
        print(f'failed to write the fasta files needed for the BlastP alignment. Exception: {e}', file=sys.stderr, flush=True)
        stop_time = time.time()
        return 2, seq_aln_dict, ids, platform.node(), worker.id, start_time, stop_time, return_code

    # run the blastp alignment
    try:
        xml_output = NcbiblastpCommandline(query=f'{aln_uuid}_query.fasta',subject=f'{aln_uuid}_sbjct.fasta',cmd=blastp_path,outfmt=5)()[0]
        blast_record = NCBIXML.read(StringIO(xml_output))
        seq_aln_dict['query_start'] = blast_record.alignments[0].hsps[0].query_start
        seq_aln_dict['query_end']   = blast_record.alignments[0].hsps[0].query_end
        seq_aln_dict['sbjct_start'] = blast_record.alignments[0].hsps[0].sbjct_start
        seq_aln_dict['sbjct_end']   = blast_record.alignments[0].hsps[0].sbjct_end
        seq_aln_dict['alignment']   = list(zip(blast_record.alignments[0].hsps[0].query,blast_record.alignments[0].hsps[0].match,blast_record.alignments[0].hsps[0].sbjct))
        seq_aln_dict['e-value']     = blast_record.alignments[0].hsps[0].expect
        seq_aln_dict['score']       = blast_record.alignments[0].hsps[0].score
        seq_aln_dict['bits']        = blast_record.alignments[0].hsps[0].bits
        seq_aln_dict['n_gaps']      = blast_record.alignments[0].hsps[0].gaps
        return_code = 0        
    except Exception as e:
        print(f'The BlastP alignment failed. Exception: {e}', file=sys.stdout, flush=True)

    # remove the fasta files created for this task
    os.remove(f'{aln_uuid}_query.fasta')
    os.remove(f'{aln_uuid}_sbjct.fasta')
    return 2, seq_aln_dict, ids, platform.node(), worker.id, start_time, time.time(), return_code


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Parsing and requesting metadata associated with structural alignment hits.')
    parser.add_argument('--input-list-file', '-inp', required=True, help='list file that contains identifiers for PDBID_CHAINID structures to be gathered.')
    parser.add_argument('--output-directory', '-outdir', required=True, help='local or global path to a directory within which structure files will be written; if not already made, the path will be made.')
    parser.add_argument('--timings-file', '-ts', required=True, help='local or global path for a CSV file for processing timings to be written; path parents need to exist already')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='local or global path for a log file to be written for this run; path parents need to exist already')
    parser.add_argument('--blastp-path', '-blastp', required=True, help='file path to the blastp executable')
    args = parser.parse_args()

    # if the user-defined output_directory is not already made, then make it
    if not Path(args.output_directory).is_dir():
        Path(args.output_directory).mkdir(mode=0o777,parents=True)

    # check to see if the output_directory has the os directory limiter
    if args.output_directory[-1] != '/':
        args.output_directory += '/'

    # start dask client.
    client = Client(name='AlignmentTaskMgr')

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger', args.output_directory + args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Timing file: {args.timings_file}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += f'Client information: {client}\n'
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')
  
    # read the input list file
    with open(args.input_list_file,'r') as list_file:
        PDBID_CHAINIDs = [line.strip() for line in list_file.readlines()]

    # parse the list to avoid redundant fetches
    pdbid_dict = {}
    for PDBID_CHAINID in PDBID_CHAINIDs:
        # assumes formatting of the strings in the list file are formatted as
        # PDBID_CHAINID...
        temp = PDBID_CHAINID.split('_')
        val = pdbid_dict.get(temp[0],[])
        val.append(temp[1])
        pdbid_dict[temp[0]] = val

    # each task is associated with one PDBID that may have one or more chain IDs
    # that will be gathered; the task list is formatted as a list of a string 
    # and a sublist. elem[0] is the PBDID accession string, elem[1] is a list of
    # chainIDs within that structure to be gathered.
    task_list = [[key,val] for key, val in pdbid_dict.items()]
    main_logger.info(f'{len(PDBID_CHAINIDs)} PDBID_CHAINID files will be gathered from {len(task_list)} unique PDB IDs.')

    # set up timing log file.
    timings_file = open(args.output_directory + args.timings_file, 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','task_type','return_code'])
    timings_csv.writeheader()
   
    # dictionary within which the structural information associated with each PDBID_CHAINID will be stored
    structure_metadata_dict = {}
    
    # list within which pdbid_chainid strings will be stored; structure files that were successfully parsed.
    pdbid_chainid_seqs = []

    # dictionary within which the pdbid_chainid key will map to the UniProt Accession ID value
    # if a pdbid_chainid does not have an associated UniprotID then a None object will be left instead
    pdbid_to_uniprot_dict= {}
    
    # list within which uniprotid strings will be stored that were successfully queried.
    uniprotid_list = []
    
    # dictionary within which the metadata associated with UniprotIDs (keys) will be stored. 
    # will be a dictionary of dictionaries, because there are many fields/types of information
    # in the Uniprot flat files. 
    uniprot_metadata_dict= {}

    # do the thing
    step0_futures = client.map(_gather_pdb_file, task_list, output_directory = args.output_directory, pure=False)
    ac = as_completed(step0_futures)
    for finished_task in ac:
        # gather the task's return objects
        step_number, results, hostname, workerid, start, stop, return_code = finished_task.result()
        # log task timing information
        append_timings(timings_csv,timings_file,hostname,workerid,start,stop,step_number,return_code)

        ### handling step 0 results:
        # if step 0 completes successfully, then save results to relevant 
        # dictionaries
        if step_number == 0 and return_code:
            # add the results dictionary to the structure_metadata_dict for 
            # storage
            structure_metadata_dict.update(results)

            # grab the structure_ids keys in the results dictionary
            structure_ids = list(results.keys())
            # once again assuming a string formatting
            pdbid = structure_ids[0].split('_')[0]
            main_logger.info(f'The structure file(s) associated with {pdbid} have been gathered, parsed, and saved to storage. Associated UniprotIDs were gathered too. nStructures: {return_code}. Took {stop-start} seconds.')

            # loop over keys and grab associated UniProtID. 
            for structure in structure_ids:
                # if no UniProtID is associated with the PDBID_CHAINID 
                # structure, or the query failed, then this will be skipped.
                uniprotid = results[structure].get('UniProtID',0)
                if uniprotid:
                    pdbid_to_uniprot_dict[structure] = uniprotid
                    # if the uniprotid has not been seen yet, create a task to 
                    # gather its metadata; otherwise, just move on
                    if uniprotid not in uniprotid_list:
                        uniprotid_list.append(uniprotid)
                        # submit a step 1 task, gathering UniProtKB flat file
                        new_future = client.submit(_query_uniprot_flat_file,uniprotid,pure=False)
                        ac.add(new_future)
        
        # if step 0 fails, write out to the log file
        elif step_number == 0:
            main_logger.info(f'The structure file(s) associated with {list(results.keys())[0]} failed to be gathered. Took {stop-start} seconds.')
            continue

        ### handling step 1 results:
        # if parsing the UniProtID flat file was successful, store the results
        # in the uniprot_metadata_dict object
        elif step_number == 1 and return_code == 0:
            uniprotid = list(results.keys())[0]
            uniprot_metadata_dict.update(results)
            main_logger.info(f'The flat file associated with {uniprotid} has been parsed. Took {stop-start} seconds.')
        # if step 1 fails, write out to the log file
        elif step_number == 1:
            uniprotid = list(results.keys())[0]
            main_logger.info(f'Requesting the flat file associated with {uniprotid} failed. Took {stop-start} seconds.')

    main_logger.info(f'Done gathering/parsing PDBID_CHAINID files and parsing UniProt flat files. Writing map and uniprot metadata dictionaries to file. {time.time()}')

    # save a structure_list file with format "path_to_structure_file nResidues"
    with open(args.output_directory + 'structures.lst','w') as out:
        for structure in structure_metadata_dict:
            file_path = structure_metadata_dict[structure].get('structure_path',0)
            last_resid = structure_metadata_dict[structure].get('struct_last_resid',0)
            first_resid = structure_metadata_dict[structure].get('struct_first_resid',0)
            nRes = last_resid - first_resid + 1
            out.write(f'{file_path} {nRes}\n')

    # save dictionary of the pdbid_chainid to uniprotid mapping
    with open(args.output_directory + 'pdbid_to_uniprotid_map.pkl', 'wb') as out:
        out_string = \
f"""{datetime.datetime.now().date()}, creator: RBD
Dictionary of mapping between PDBID_CHAINID and UniProtKB accession IDs. Format:
    0) string, this statement,
    1) dict, keys are PDBID_CHAINID strings and values are UniProtKB accession IDs."""
        pickle.dump([out_string,pdbid_to_uniprot_dict],out)

    # save dictionary of uniprot accession id meta data
    with open(args.output_directory + 'uniprot_metadata.pkl', 'wb') as out:
        out_string = \
                f"""{datetime.datetime.now().date()}, creator: RBD
                Stash of UniProtKB flat file metadata. Format:
                    0) string, this statement,
                    1) dict, keys are UniProtKB accession ID strings and values are subdictionaries filled with metadata information gathered from the UniProtKB ID flat file.\n\t""" + "Ex: key 'S4UX02' has value {'entry_name': 'CYPH1_SALMI','status': 'Reviewed','sequence_length': 495, ... }\n\n"+"""
                
                Subdictionary organization:
                    'entry_name': qualitative name for the UniProtKB entry.
                    'status': qualitative description of how well studied/reviewed the entry is.
                    'sequence_length': number of residues in the sequence reported in 'sequence'.
                    'original_date': date for first entry made public (format: YYYY-MM-DD).
                    'sequence_date': date for sequence being made public or changed (format: YYYY-MM-DD).
                    'modified_date': date for last change made to entry information (format: YYYY-MM-DD).
                    'descriptions': list of strings, each string is a 'DE' line in the UniProtKB flat file.
                    'gene_name': common name for the encoding gene associated with the entry's protein.
                    'species_name': originating species' name for the entry's protein.\n\t'organelle_loc': cell localization information.
                    'taxon': taxonomy information.
                    'taxon_id': taxonomy identification string.
                    'database_references': list of lists, where each sublist is the metadata gathered for one of the UniProtKB entry's categorization within a database; many databases are considered; see formatting of 'DR' lines for more details.
                    'primary_evidence': qualitative string that describes to what level the experimental evidence reports on this entry.
                    'keywords': strings that describe this entry.
                    'features': list of lists, where each sublist is the metadata gathered for one residue-level piece of information about the UniProtKB entry; many different types of feature types; see formatting of 'FT' lines. for more details.
                    'sequence': amino acid sequence string, including the fasta-like header line.
                    'ecIDs': list of strings, empty if no EC ID(s) are found in the UniProtKB flat file.
                    'query_status': return code from the python request call for gathering the UniProtKB flat file."""
        pickle.dump([out_string,uniprot_metadata_dict],out)

    main_logger.info(f'Beginning to run blastp alignments between structure and uniprot sequences to get residue index mapping. {time.time()}')
    
    
    task2_futures = []
    # only loop over pdbid_chainids that were successfully parsed by step 0
    for pdbid_chainid, uniprotid in pdbid_to_uniprot_dict.items():
        # check to see that the pdbid_chainid maps to a real uniprotid
        if uniprotid:
            main_logger.info(f'{pdbid_chainid} maps to {uniprotid}. Submitting a blastp task to align the two sequences.')
            # submit a task to run the _blast_aln function, aligning structure
            # sequence to the respective UniProt flat file's sequence
            struct_seq = structure_metadata_dict[pdbid_chainid]['struct_seq']
            uniprot_seq = uniprot_metadata_dict[uniprotid]['sequence']
            new_future = client.submit(_blast_aln, struct_seq, uniprot_seq, args.blastp_path, [pdbid_chainid, uniprotid], pure=False)
            task2_futures.append(new_future)

    task2_ac = as_completed(task2_futures)
    for finished_task in task2_ac:
        task_num, results, ids, hostname, workerid, start, stop, return_code = finished_task.result()
        pdbid_chainid = ids[0]
        uniprotid     = ids[1]
        
        # log task timing information
        append_timings(timings_csv,timings_file,hostname,workerid,start,stop,task_num,return_code)
        # if the sequence alignment was successful
        if return_code == 0:
            main_logger.info(f'The blastp alignment run between {pdbid_chainid} and {uniprotid} sequences has completed. Return code: {return_code}. Took {stop - start} seconds.')
            
            structure_metadata_dict[pdbid_chainid]['UniProt_blast_aln'] = results
            structure_metadata_dict[pdbid_chainid]['UniProt_features'] = []
            ### logic:
            # features are lists of 3+ elements; element 1 holds the indices of
            # residue(s) that describe said feature; these resids are in seq space
            # of the uniprot sequence and are 1-indexed.
            # Need to check that a feature is within the aligned region of the two
            # sequences by boolean testing if feat_resid(s) are in 
            # range(sbjct_start,sbjct_end+1); 
            # all indices at the moment are 1-indexed so no need to worry about
            # shifting indices.
            # seq_aln_delta: mapping shift for sequence alignment; 
            #                sbjct_start to query_start 
            # struct_delta:  mapping shift from struct_seq to structure_resids
            #                struct_seq was pulled from the pdb file so first index
            #                is 1 and maps to the residues resid in the pdb file. 
            seq_aln_delta = results['query_start'] - results['sbjct_start']
            struct_delta  = structure_metadata_dict[pdbid_chainid]['struct_first_resid'] - 1
            for feature in uniprot_metadata_dict[uniprotid]['features']:
                # expected format is either e.g. "49" or "49..63"
                try:
                    flat_file_resindex = [int(index) for index in feature[1].split('..')]
                # sometimes you might see e.g. "?..49" or "<49", we are gonna 
                # ignore these features
                except:
                    main_logger.info(f'{uniprotid}: {feature} was not parsed.')
                    continue
                # flat_file_resindex can be a list of one or two ints; if any of 
                # these ints is within the range(sbjct_start,sbjct_end), then we are
                # interested in what that feature is and wanna convert its residue
                # range to a structure residue indexed range
                # this conversion requires a couple hops to do:
                # the aln's sbjct resid to the aln's query resids (all 1-indexed)
                # and then
                # from aln's query resids to the structure's residue-indexing, which
                # is not necessarily one-indexed
                if sum([seqindex in list(range(results['sbjct_start'],results['sbjct_end']+1)) for seqindex in flat_file_resindex]) >= 1:
                    struct_resindex = [resindex + seq_aln_delta + struct_delta for resindex in flat_file_resindex]
                    temp_feature = copy.deepcopy(feature)
                    # these struct_resindex numbers may not correspond well to the 
                    # structure_seq because that sequence is naive of missing 
                    # residues (whether unresolved in the structure or 
                    # non-standard); these missing residues will appear as gaps in
                    # the query alignment elements from the list of tuples
                    temp_feature[1] = struct_resindex
                    structure_metadata_dict[pdbid_chainid]['UniProt_features'].append(temp_feature)
        else:
            main_logger.info(f'The blastp alignment run between {pdbid_chainid} and {uniprotid} sequences failed. Return code: {return_code}. Took {stop - start} seconds.')

    # save dictionary of pdbid_Chainid metadata dictionary
    with open(args.output_directory + 'pdbid_chainid_metadata.pkl', 'wb') as out:
        output_string = f"""
        {datetime.datetime.now().date()}, creator: RBD
        Dictionary where keys are PDBID_CHAINID strings and values are a dictionary filled with information about the structure.
        Keys and Values in subdictionaries:
            'structure_path': string to the path for the .pdb file for this structure; update if structures are moved.
            'struct_first_resid': first residue number in the structure.
            'struct_last_resid': last residue number in the structure.
            'struct_seq': amino acid sequence for resolved structure residues.
            'UniProtID': string, if the structure maps to a UniProt Accession ID, then the value is that ID.
            'UniProt_features': list of lists, each list is a residue-level feature associated with the UniProt entry that has now been mapped to the structure. See below for description of data organization.
            'UniProt_blast_aln': dictionary with keys being descriptive strings for blast results, values being those results.

        UniProt_features list description:
        ['TypeOfFeature',[FirstResidueNumber,LastResidueNumber],'Further','Descriptive','Strings']
            The 0th element in each list is the type of feature.
            The 1st element in each list is a list of one or two integers; if only one, then the feature is only associated with that residue; if two integers, then the feature is associated with all residues in the inclusive range of residue numbering.
            Further list elements are descriptive strings for the feature.
        Example:['ACT_SITE',[81],'/note="Schiff-base intermediate with D-ribose 5-phosphate','/evidence=ECO:0000255|HAMAP-Rule:MF_01824,','ECO:0000305|PubMed:15911615"']

        UniProt_blast_aln dictionary description:
            'query_start', 'query_end': first and last residue index associated with the sequence of residues listed in the structure file.
            'sbjct_start', 'sbjct_end': first and last residue index associated with the UniProt Accession sequence that was aligned during the BlastP alignment.
            'alignment': list of tuples, where each tuple is the alignment residue pairs output from BlastP.
            'e-value': e-value output from BlastP; should be extremely small since we are aligning equivalent sequences.
            'score': score value output from BlastP; should be large since we are aligning equivalent sequences.
            'bits': bits value output from BlastP; should be large since we are aligning equivalent sequences.
            'n_gaps': number of gaps used to optimize sequence alignment; should be small but may have gaps if there are unresolved residues in the structure."""
        pickle.dump([output_string,structure_metadata_dict],out)

    # close log files and shut down the cluster.
    timings_file.close()
    client.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

