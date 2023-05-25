#!/usr/bin/env python3
"""
    Query RCSB site Uniprot ID associated with a given protein and chain ID.

    (Written by Mark Coletti, colettima@ornl.gov.)
"""
import requests

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
