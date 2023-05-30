#!/usr/bin/env python3
"""
    Request Uniprot flat files associated with given UniProt Accession IDs

    (Written by Russell B Davidson, davidsonrb@ornl.gov, rbdavid)
"""

import requests
import datetime
import copy
import re

flat_file_url = 'https://www.uniprot.org/uniprot/%s.txt'
dictionary_namespace = {'entry_name': '',           # ID line
                        'status': '',               # ID line
                        'sequence_length': '',      # ID line
                        'original_date': '',        # 1st DT line
                        'sequence_date': '',        # 2nd DT line
                        'modified_date': '',        # 3rd DT line
                        'descriptions': [],         # DE lines
                        'gene_name': [],            # GN lines
                        'species_name': '',         # OS line
                        'organelle_loc':'',         # OG line
                        'taxon':'',                 # OC line
                        'taxon_id':'',              # OX line
                        'database_references': [],  # DR lines
                        'primary_evidence': '',     # PE line
                        'keywords': [],             # KW lines
                        'features': [],             # FT lines
                        'sequence': '>',            # SQ and blank lines afterwards
                        'ecIDs': [],                # Gather any lines that have EC IDs in them
                        'query_status': 0}

def request_uniprot_metadata(accession_id):
    """
    """
    try:
        response = requests.get(flat_file_url %(accession_id))
    except Exception as e:
        print(f'Requesting the associated flat file for {accession_id} failed.')
        return {'query_status':-9999}
    
    status_code   = response.status_code
    response_text = response.text
   
    # successful response
    if status_code == 200:
        uni_dict = copy.deepcopy(dictionary_namespace)
        uni_dict['query_status'] = 200
        response_lines = response_text.split('\n')
        for line in response_lines:
            # end of file is denoted by '//' so break from the for loop if it occurs
            if line == '//':    # used to denote end of file for the UniProt flat files
                break
            # parse first line, the identity line
            # multiple elements per line
            elif line[:2] == 'ID':
                temp = re.findall('(\w+)',line)
                uni_dict['entry_name'] = temp[1]
                uni_dict['status'] = temp[2]
                uni_dict['sequence_length'] = int(temp[3])
            # parse the chronology lines; each entry has three DT lines
            # one date per line
            elif line[:2] == 'DT':
                dt = re.findall('([\w-]+)',line)[1]
                if 'integrated' in line:
                    uni_dict['original_date'] = datetime.datetime.strptime(dt,'%d-%b-%Y').strftime('%Y-%m-%d')
                elif 'sequence version' in line:
                    uni_dict['sequence_date'] = datetime.datetime.strptime(dt,'%d-%b-%Y').strftime('%Y-%m-%d')
                elif 'entry version' in line:
                    uni_dict['modified_date'] = datetime.datetime.strptime(dt,'%d-%b-%Y').strftime('%Y-%m-%d')
            # parse the description lines; gathering all DE lines
            elif line[:2] == 'DE':
                #description = re.findall('(?<=Full\=)[^,;]*|(?<=Short\=)[^,;]*|(?<=EC\=)[^,;]*',line)
                uni_dict['descriptions'].append(line[2:].strip())
                # gathering DE line instances where  EC IDs are reported
                if 'EC=' in line:
                    uni_dict['ecIDs'] += re.findall('(?<=EC\=)[^,;\s]*',line)
            # parse the gene name line(s)
            # multiple names or synonyms possible per line
            elif line[:2] == 'GN':
                names = re.findall('(?<=[\=,])[^,;]*',line)
                for name in names:
                    uni_dict['gene_name'].append(name)
            # parse the species name line(s); just add the names to a string, no need to parse further
            elif line[:2] == 'OS':
                uni_dict['species_name'] += line[5:].strip()
            # parse the organelle line to collect localization info
            elif line[:2] == 'OG':
                uni_dict['organelle_loc'] += line[5:].strip()
            # parse the database reference lines
            # these DR lines are extremely diverse in format and content... 
            # got lazy with regex and just parsing into a list of lists, just 
            # separating by the semicolon delimiter used for these lines;
            # 0th element of each sublist will be the database name
            # further elements are accession ids or other metadata
            elif line[:2] == 'DR':
                dr = [elem.strip() for elem in re.findall('[^;\n]+',line[5:].rstrip('.'))]
                uni_dict['database_references'].append(dr)
            # parse the primary evidence line
            # just get the number, always 6th character in line
            elif line[:2] == 'PE':
                uni_dict['primary_evidence'] = line[5]
            # parse the KW lines
            # following similar lazy behavior as DR lines since keywords can be 
            # very diverse in format/length
            elif line[:2] == 'KW':
                kw = [elem.strip() for elem in re.findall('[^;\n]+',line[5:].rstrip('.'))]

            # parse the features lines
            # these FT lines are extremely diverse in format and content...
            # if the line is associated with a new feature, the 5th character
            # will be the first letter of the feature type and not a ' '
            # in this case, the last element of the split line will be an 
            # integer or a pair of integers; these are the first and last
            # sequence indices associated with the feature; parse these numbers
            # append [feature_type, residue_range] to the features dict list
            elif line[:2] == 'FT' and line[5] != ' ':
                uni_dict['features'].append(line[5:].split())
            # if the line's 5th character does equal ' ', then this line is 
            # associated with a previous feature
            # in this case, append the line's information to the previous list.
            # these lines are usually notes, ids, and evidence information.
            # just strip new line characters and append the string to the end 
            # of the list
            elif line[:2] == 'FT' and line[5] == ' ':
                uni_dict['features'][-1].append(line[2:].strip())
            # parse the sequence lines
            # the 'SQ' key is the notification that the sequence information
            # begins on this line
            elif line[:2] == 'SQ' and line[5] != ' ':
                uni_dict['sequence'] += line[5:] + '\n'
            # the only instance where the first two characters of a line == '  '
            # is a line that continues with sequence information
            # of course, these lines are formatted; blocks of 10 AAs separated
            # by a space character; 60 AAs on a line
            elif line[:2] == '  ':
                temp = line.split()
                seq_string = ''.join(temp)
                uni_dict['sequence'] += seq_string
        
        return uni_dict

    # if the request response fails with known 'failure' status codes
    elif status_code in [400,404,410,500,503]:
        print(uniprotID, status_code, response_text)
        return {'query_status': status_code}
    # if the request response fails for any other reason
    else:
        print(uniprotID, status_code, 'Something really funky is going on')
        return {'query_status': status_code}


if __name__ == '__main__':
    # real uniprot accession id
    results = request_uniprot_metadata('A0A0M3KL33')
    assert results['status'] == 0

    # fake uniprot accession id
    results = request_uniprot_metadata('BOGUS')
    assert results['status'] != 0

