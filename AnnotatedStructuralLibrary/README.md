# Structure Library and Metadata Gathering
*Written by Russell B. Davidson and Mark Coletti*
Version 1.0.0

Bash and python scripts to automate the download of structure files from the 
RCSB database and gather the associated metadata from the UniProtKB database 
when available. Residue-level feature metadata is transferred to the associated
structure. This metadata is then used to garner insights from structural 
alignments against the annotated library of structures developed from this code.

## Using the Workflow



### Workflow Steps:
1) Download the structural data associated with PDB IDs, separate the whole 
structural data into individual chains, remove non-standard residues, and save 
to storage.
2) Parse the structure files to gather sequences and residue indices.
3) Query the RCSB through their GraphQL API to gather the UniProt Accession IDs, 
when available. 
4) Query the UniProtKB API to gather and parse the metadata flat file associated 
with each Accession ID. 
5) Run a BlastP sequence alignment analysis to map the UniProtKB sequence to the 
structure sequence. Once this map is made, feature elements from the UniProtKB 
flat files can be transferred to the structure's residues. 


## BlastP Installation Instructions
Instructions for installing the BlastP command line application are available 
through the [NCBI Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK569861/).


## Container Installation 



