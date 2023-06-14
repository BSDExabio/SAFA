# StructuralAlignmentWorkflow
*Written by Russell B. Davidson and Mark Coletti*, Version 1.1.0

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
See the `example` subdirectory that houses two instances of the structural 
alignment workflow where the total number alignments are 'small' and 'large'.
The two bash scripts were used to run the workflow for these two examples. 

## Installation of US-align2

US-align2 is easily installed via the associated github repository. 

```
git clone https://github.com/pylelab/USalign.git
```

That repository can also be updated to the most current version by pulling from 
the remote repository: 

```
git pull
```

## Parsing US-align2 output
Functions used to parse the US-align2 results are provided in the 
`USalignParser` subdirectory. A jupyter notebook is provided here to demonstrate
how that set of functions are used in the Alignment workflow to gather the 
relevant quantitative results for each alignment calculation. 


## Post-processing and analysis
Scripts are provided in the `PostProcessingAnalyses` subdirectory that are used
to take an ensemble of alignment results and provide novel insights, such as 
enzyme commission (EC) number annotation hypotheses or visualization codes for 
high quality alignment results. These codes are often dependent on the 
structural library with which the structural alignment workflow utilized. See 
the PDB70 structural library 
([https://doi.org/10.5281/zenodo.7953087](https://doi.org/10.5281/zenodo.7953087))
data set for the structural library used to annotate enzymes in the *S. divinum*
structural proteome. 
