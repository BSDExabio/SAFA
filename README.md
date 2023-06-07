# SAFA - Structural Alignment for Functional Annotation

This repository contains the scripts and workflow codes used to perform structural alignment of protein model(s) against a library of structures (e.g., the PDB70 structural library), along with parsers and post-processing/analysis scripts (StructuralAlignmentWorkflow). Additionally, the processing codes used to gather the metadata associated with a structural library are included to enable users to develop their own annotated structural libraries (AnnotatedStructuralLibrary). For both workflows, jupyter notebooks are included in addition to the codes to enable demonstration the steps of each workflow and visualization of results. Tests for each codebase are provided as is a small scale example of all components of SAFA. 

# Supplementary Datasets
The codes provided here are supplemented by the PDB70 structural library ([https://doi.org/10.5281/zenodo.7953087](https://doi.org/10.5281/zenodo.7953087)), the ModelArchive structural repository for *S. divinum* proteome ([https://modelarchive.org/doi/10.5452/ma-ornl-sphdiv](https://modelarchive.org/doi/10.5452/ma-ornl-sphdiv)), and the functional annotation dataset for the *S. divinum* proteome ([https://doi.org/10.13139/ORNLNCCS/1973885](https://doi.org/10.13139/ORNLNCCS/1973885)). 

1. The PDB70 structural library represents a set of experimentally-determined protein structures that were used as the structural alignment queries for the functional annotation of the *S. divinum* proteome. The Zenodo dataset also includes the metadata for these strucutures, gathered from the UniProtKB database. 

2. For each protein in the *S. divinum* proteome, the highest pTM scoring predicted structural model is deposited in the ModelArchive repository. The pLDDT and PAE metrics for these top-ranked structural models are also available. 

3. Finally, the sequence and structural alignment results are available through Constellation, the Oak Ridge Leadership Computing Facility (OLCF) DOI-based science network for supercomputing data; this data repository also houses all five structural models output from AlphaFold. 

