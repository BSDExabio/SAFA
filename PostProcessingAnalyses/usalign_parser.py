#!/usr/bin/env python3
"""
    For parsing USalign alignment files (-outfmt 0 or -1)

"""
import time
from pathlib import Path
import pandas as pd
import re


def parse_usalign_file(fn,alignment_type):
    """ Read given USalign alignment file

    :param fn: path string for alignment log file
    :param alignment_type: string parameter to set format of alignment output to parse; currently accepted: CP, SNS, or FNS
    :return: dictionary of quantitative results associated with the alignment. Keys:
                    'struct1': string, path to the query structure; always assumed to be the model structure
                    'struct2': string, path to the target structure; always assumed to be the template structure
                    'struct1_chainID': string, chain ID string associated with the query structure
                    'struct2_chainID': string, chain ID string associated with the target structure
                    'tmscore1': float, TMscore value normalized by the query structure's length
                    'tmscore2': float, TMscore value normalized by the target structure's length
                    'Len1': float, the query structure's length
                    'Len2': float, the target structure's length
                    'd0_1': float, the normalized distance for query structure; used in TMscore metric calculation
                    'd0_2': float, the normalized distance for target structure; used in TMscore metric calculation
                    'map_1_to_2': dictionary, keys are query structure's ('1') residue indices (string) that map to a tuple with query residue name, target structure's ('2') residue index (string), target residue name, and alignment distance (float, only available for SNS and FNS alignments); optional, will not be present if alignment_type != 'CP', 'SNS', 'FNS'
    """
    start_time = time.time()
    results = {}
    with open(fn,'r') as aln_log_file:
        lines = aln_log_file.readlines()

    ### gather relevant lines common across all USalign outputs
    # query structure lines
    struct1_lines = [line.strip() for line in lines if 'Structure_1:' in line]
    results['struct1'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct1_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    results['struct1_chainID'] = re.search(r'(?<=[:])\w+',struct1_lines[0])[0]    # finds the chainID that was used in the alignment
    results['tmscore1'], results['Len1'], results['d0_1'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct1_lines[2])]    # gathering quantitative metrics associated with TM-score, Length of query structure, and d0 value for the alignment
    
    # target structure lines
    struct2_lines = [line.strip() for line in lines if 'Structure_2:' in line]
    results['struct2'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct2_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    results['struct2_chainID'] = re.search(r'(?<=[:])\w+',struct2_lines[0])[0]    # finds the chainID that was used in the alignment
    results['tmscore2'], results['Len2'], results['d0_2'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct2_lines[2])]    # gathering quantitative metrics associated with TM-score, Length of target structure, and d0 value for the alignment
    
    # Alignment overview line
    aln_lines    = [line.strip() for line in lines if 'Aligned length=' in line]
    results['alnLen'], results['rmsd'], results['alnSeqID'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',aln_lines[0])]    # gathering quantitative metrics associated with Length of Alignment, RMSD, and sequence ID for the aligned residues

    if alignment_type.upper() == 'CP':
        # gather the alignment mapping
        map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]

        # dictionary of tuples; keys are mobile residue index with values being a tuple of (target residue index, mobile resname, target resname)
        results['map_1_to_2'] = {}
        for line in map_lines:
            # line format: 'CA  LEU A 111 \t CA  GLN A  97'
            # awkward white spaces:      ^^
            # wanna keep: {'111': ('97','LEU','GLN')}
            temp = [elem.strip() for elem in line.split('\t')]
            results['map_1_to_2'].update({temp[0][-4:].strip(): (temp[1][-4:].strip(),temp[0][4:7],temp[1][4:7])})
    
    elif alignment_type.upper() in ['SNS','FNS']:
        # gather the alignment mapping
        map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]
        
        # dictionary of tuples; keys are mobile residue index with values being a tuple of (target residue index, mobile resname, target resname, distance)
        results['map_1_to_2'] = {}
        for line in map_lines:
            # line format: 'CA  LEU A 111 \t CA  GLN A  97 \t    1.568'
            # awkward white spaces:      ^^               ^^
            # wanna keep: {'111': ('97','LEU','GLN',1.568)}
            temp = [elem.strip() for elem in line.split('\t')]
            results['map_1_to_2'].update({temp[0][-4:].strip(): (temp[1][-4:].strip(),temp[0][4:7],temp[1][4:7],float(temp[2]))})

    else:
        print(f"'{alignment_type}' not expected by parser function. Only returning alignment quantitative metrics. No mapping.")
        stop_time = time.time()
        return results, start_time, stop_time, 1

    stop_time = time.time()
    return results, start_time, stop_time, 0


if __name__ == '__main__':
    # test harness for parsers

    result_dict, start, stop, return_code = parse_tmalign_score_file('/home/russ/Projects/ornl_DLFA/structural_alignment_codes/USalign_results/3cna_2pel_usalign_tmalign.log','None')

    result_dict, start, stop, return_code = parse_tmalign_score_file('/home/russ/Projects/ornl_DLFA/structural_alignment_codes/USalign_results/3cna_2pel_usalign_cp.log','CP')

    result_dict, start, stop, return_code = parse_tmalign_score_file('/home/russ/Projects/ornl_DLFA/structural_alignment_codes/USalign_results/3cna_2pel_usalign_sNS.log','sNS')

    result_dict, start, stop, return_code = parse_tmalign_score_file('/home/russ/Projects/ornl_DLFA/structural_alignment_codes/USalign_results/3cna_2pel_usalign_fNS.log','fNS')

    pass
