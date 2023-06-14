# Structure Library and Metadata Gathering
*Written by Russell B. Davidson and Mark Coletti*
Version 1.1.0

Bash and python scripts to automate the download of structure files from the 
RCSB database and gather the associated metadata from the UniProtKB database, 
when available. Residue-level feature metadata is transferred to the associated 
structure. This metadata can then be used in a broad range of analyses such as 
functional annotation.  

## Using the Workflow

The workflow script (`processing_pdb_library.py`) has five input variables that
must be defined:

1) --input-list-file, -inp: a list of RCSB PDB identifiers as well as the chain 
identifier. Each line should have the format of `XXXX_Y` where XXXX is the PDB 
ID and Y is the chain identifier. 

2) --output-directory, -outdir: a local or global path pointing towards a 
directory within which all files will be written. 

3) --timings-file, -ts: a file name to be used for writing the timing 
information for each task in the dask.distributed workflow. 

4) --tskmgr-log-file, -log: a file name to be used for writing the loggging
messages for the workflow. 

5) --blastp-path, -blastp: a local or global path pointing to the BlastP 
executable. 


The only user-prepared input necessary for running this workflow is a list of 
RCSB PDB identifiers as well as the chain identifier. An example of such a file 
is provided in the `example` subdirectory.


### Workflow Steps:
1) Download the structural data associated with PDB IDs, separate the whole 
structural data into individual chains, remove non-standard residues, and save 
to storage. Parse the structure files to gather sequences and residue indices.
Query the RCSB through their GraphQL API to gather the UniProt Accession IDs, 
when available. 
2) Query the UniProtKB API to gather and parse the metadata flat file associated 
with each Accession ID. 
3) Run a BlastP sequence alignment analysis to map the UniProtKB sequence to the 
structure sequence. Once this map is made, feature elements from the UniProtKB 
flat files can be transferred to the structure's residues. 

### Simple Example
See the `example` subdirectory that contains the input list as well as all 
output from the workflow associated with that list. The 
`processing_pdb_library.sh` is the bash script that was used to run this 
example.

## BlastP Installation Instructions
Instructions for installing the BlastP command line application are available 
through the [NCBI Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK569861/).


