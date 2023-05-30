# StructuralAlignmentWorkflow
*Written by Russell B. Davidson and Mark Coletti*
Version 1.0.0

A dask.distributed workflow code to automate the alignment of a library of 
structures against a library of structures. This code was used to automate the 
structural alignment for functional annotation reported in (CITE). Currently,
the workflow utilizes US-align2 (https://github.com/pylelab/USalign) to
perform the single-chain to single-chain protein structural alignment 
calculations. 

## Using the Workflow

### Workflow Steps:
1) Gather the lists of structures to be aligned. A set of query structures (ex: 
the PDB70 structural library) and a set of target structures (ex: the structral 
models of the *S. divinum* proteome) are input files to the workflow. These 
lists are sorted largest structures to smallest based on the number of residues 
in the structure. This organization is performed as a quick and dirty sorting 
of tasks that has been seen to boost efficiency.  
2) The dask.distributed workflow submits the alignment tasks to the worker 
resources. If the number of tasks is small, this is brute forced where each 
alignment between query and target structures are a single task. Alternatively, 
if the number of alignments to be performed is large, then the query list is 
partitioned into sublists. Each task then represents alignments between the 
target model and a list of query structures.
3) Results of all alignments are gathered in temporary output files and parsed 
for ranking. During this process, file IO of strong alignments is performed if 
the user defined the metric and cutoff for defining what an alignment hit was. 
The full set of results are written to a target-specific output directory.  

### Simple Example


## Container Installation


## Installation of US-align2



