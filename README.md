# SAFA - Structural Alignment for Functional Annotation

This repository contains the scripts and workflow codes used to perform structural alignment of protein model(s) against a library of structures (e.g., the PDB70 structural library), along with parsers and post-processing/analysis scripts (StructuralAlignmentWorkflow). Additionally, the processing codes used to gather the metadata associated with a structural library are included to enable users to develop their own annotated structural libraries (AnnotatedStructuralLibrary). 
<!---
For both workflows, jupyter notebooks are included in addition to the codes to enable demonstration the steps of each workflow and visualization of results. 
-->
Tests for each codebase are provided as is a small scale example of all components of SAFA. 

# DATASETS for *Sphagnum Divinum* paper
For each protein in the *Sphagnum divinum* proteome, the highest pTM scoring predicted structural model is deposited in the ModelArchive repository at https://modelarchive.org/doi/10.5452/ma-ornl-sphdiv. 

The sequence and structural alignment results are available through Constellation, the Oak Ridge Leadership Computing (OLCF) DOI-based science network for supercomputing data, at (https://doi.org/10.13139/ORNLNCCS/1973885); this data repository also houses all five structural models output from AlphaFold. 

The PDB70 structural library and associated metadata is available at https://doi.org/10.5281/zenodo.7953087. 

