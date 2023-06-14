#!/usr/bin/env python3
"""
    Query RCSB through various APIs to gather structures and metadata such as 
    Uniprot ID associated with a given protein and chain ID.

    (Written by Russell Davidson and Mark Coletti)
"""
import sys
import datetime
import requests
import mmtf
import warnings
import MDAnalysis

QUERY_STR = \
'''
query($id: String!)
{
  entries(entry_ids:[$id]){
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
      }
      uniprots {
        rcsb_id
      }
    }
  }
}
'''


def gather_structure(PDBID):
    """
    There are multiple API/access strategies that we can use; here, the mmtf 
    module is used to gather the structure. If this fails, no structure object 
    is returned.

    NOTE: Implement other modules or API strategies for downloading structures
    and associated data.
    
    :param PDBID: string, RCSB accession ID for a specific structure.
    :return: a MDAnalysis Universe object associated with the requested 
             structure.
    """
    try:
        # fetch the whole PDBID entry using mmtf
        pdb_obj = mmtf.fetch(PDBID)
        # gather relevant meta info about the structure
        PRODUCER  = f'RBD {datetime.datetime.now().date()}'
        STRUCTURE = pdb_obj.structure_id
        VERSION   = pdb_obj.mmtf_version
        TITLE     = pdb_obj.title
        DEPDATE   = pdb_obj.deposition_date
        RELDATE   = pdb_obj.release_date
        EXPMETHOD = pdb_obj.experimental_methods
        SPCGRP    = pdb_obj.space_group
        R_FREE    = pdb_obj.r_free
        R_WORK    = pdb_obj.r_work
        RESOLUTION= pdb_obj.resolution
    except Exception as e:
        print(PDBID, 'failed to grab mmtf obj', str(e), file=sys.stderr, flush=True)
        return None

    try:
        with warnings.catch_warnings():
            # ignore some annoying warnings from sel.write line due to missing information (chainIDs, elements, and record_types). 
            warnings.simplefilter('ignore',UserWarning)

            # load the mmtf object into MDAnalysis
            u = MDAnalysis.Universe(pdb_obj)
            
            # need to check for non-standard amino acid resnames that MDAnalysis 
            # cannot parse into a sequence.
            # replace these nonstandard resnames with their closest analog 
            # ASN1 -> ASN, CYM  -> CYS, HYP -> PRO, MSE -> MET, PYL -> LYS, 
            # SEC  -> CYS; PGLU -> GLU 
            # list may be incomplete and so this may fail due to strange cases
            # standard "protein" residues that MDAnalysis comprehends: 
            # https://userguide.mdanalysis.org/1.1.1/standard_selections.html
            # instances where this code will fail: 
            # - C or N termini residues with the termini as a prefix to the resname
            sel = u.select_atoms('protein and not resname ACE ALAD ASF CME CYX DAB NME QLN ORN')
            for resid in range(sel.n_residues):
                if sel.residues[resid].resname == 'ASN1':
                    sel.residues[resid].resname = 'ASN'
                elif sel.residues[resid].resname == 'CYM':
                    sel.residues[resid].resname = 'CIS'
                elif sel.residues[resid].resname == 'HIP':
                    sel.residues[resid].resname = 'HIS'
                elif sel.residues[resid].resname == 'HISD':
                    sel.residues[resid].resname = 'HIS'
                elif sel.residues[resid].resname == 'HISE':
                    sel.residues[resid].resname = 'HIS'
                elif sel.residues[resid].resname == 'HSP':
                    sel.residues[resid].resname = 'HIS'
                elif sel.residues[resid].resname == 'HYP':
                    sel.residues[resid].resname = 'PRO'
                elif sel.residues[resid].resname == 'MSE':
                    sel.residues[resid].resname = 'MET'
                elif sel.residues[resid].resname == 'PGLU':
                    sel.residues[resid].resname = 'GLU'
                elif sel.residues[resid].resname == 'PYL':
                    sel.residues[resid].resname = 'LYS'
                elif sel.residues[resid].resname == 'SEC':
                    sel.residues[resid].resname = 'CYS'

            # save meta information as remarks lines that will be written out to PDB
            u.trajectory.remarks = [f'FILE PRODUCER: {PRODUCER}',
                                    f'MMTF VERSION: {VERSION}',
                                    f'STRUCTURE TITLE: {TITLE}',
                                    f'STRUCTURE PDBID: {STRUCTURE}',
                                    f'DEPOSITION DATE: {DEPDATE}',
                                    f'RELEASE DATE: {RELDATE}',
                                    f'EXPERIMENTAL METHOD: {EXPMETHOD}',
                                    f'SPACE GROUP: {SPCGRP}',
                                    f'R_FREE: {R_FREE}',
                                    f'R_WORK: {R_WORK}',
                                    f'RESOLUTION: {RESOLUTION}']
            return u
    except Exception as e:
        print(PDBID, 'failed to load/alter the pdb_obj into MDAnalysis', str(e), file=sys.stderr, flush=True)
        return None


def query_uniprot_str(protein_chain):
    """ Get Uniprot ID for the given protein and protein chain

    Just a convenience wrapper for query_uniprot().

    :param protein_chain: String that contains "AAAA_C"
    :return: Uniprot ID or None if not found
    """
    protein, chain = protein_chain.split('_')
    return query_uniprot(protein, chain)


def query_uniprot(protein, chain):
    """ Get Uniprot ID for the given protein and protein chain

    This will look for the chain ID in the current and any past chain IDs.

    :param protein: for which to get Uniprot ID
    :param chain: specific chain in that protein
    :return: Uniprot ID or None if not found
    """
    uniprot_id = None

    r = requests.post('https://data.rcsb.org/graphql',
                      json={'query' : QUERY_STR,
                            'variables' : {'id' : protein.upper()}})
    data = r.json()

    if 'error' in data:
        # There was something wrong with the query
        # TODO extract some meaningful error from this for an informative
        #  exception
        return None

    ### NOTE: THIS LINE IS THE CAUSE OF A FEW BUGS IN THE ORIGINAL PDB70 DATASET
    chain = chain.upper() # Normalize to chain IDs being upper case

    for entry in data['data']['entries']:
        for entity in entry['polymer_entities']:
            if chain in entity['rcsb_polymer_entity_container_identifiers']['asym_ids'] or chain in entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']:
                if entity['uniprots'] != None:
                    uniprot_id = entity['uniprots'][0]['rcsb_id']
                    break
                else:
                    continue

    return uniprot_id


if __name__ == '__main__':
    # Simple single chain
    results = query_uniprot('6lzm', 'A')
    assert results == 'P00720'

    # Test of convenience function
    results = query_uniprot_str('6lzm_A')
    assert results == 'P00720'

    # More than one chain
    results = query_uniprot('5fvk', 'B')
    assert results == 'P52917'

    # Finding among many chains
    results = query_uniprot('4v8m', 'C')
    assert results == 'Q385D9'

    # Finding using old chain name
    results = query_uniprot('4v8m', 'A2')
    assert results == 'Q385D9'

    # Check for failed search for chain that does not exist
    results = query_uniprot('6lzm', 'ZZ')
    assert results == None

    # Check for entity that has no uniprot ID associated with it
    results = query_uniprot('5t0w', 'C')
    assert results == None

    # Check for failed search for protein that does not exist
    results = query_uniprot('BOGUS', 'A')
    assert results == None


    print('All tests passed.')
