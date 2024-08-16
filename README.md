# The Germline Genomics of Cancer (G2C)
## Variant calling from germline whole-genome sequencing (WGS) across cancer types

Copyright (c) 2023-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu), [Noah Fields](mailto:noah_fields@dfci.harvard.edu), and the Van Allen, Gusev, and Haigis laboratories at Dana-Farber Cancer Institute.  
Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  

#### _Note: this repository is under active development. More documentation will be added as the project evolves._

---  

## Synopsis    

This repository contains the working code and scripts used to detect, genotype, filter, and annotate germline variants from germline WGS across cancer types  

---  

## Table of Contents  

| Directory | Description |  
| :--- | :--- |  
| [`docker/`](https://github.com/talkowski-lab/dsmap/tree/main/docker) | Instructions for building project-related Docker images |   
| [`scripts/`](https://github.com/talkowski-lab/dsmap/tree/main/scripts) | Stand-alone scripts called by various workflows |   
| [`shell/`](https://github.com/talkowski-lab/dsmap/tree/main/shell) | Shell snippets for running specific processes or analyses |   
| [`wdl/`](https://github.com/talkowski-lab/dsmap/tree/main/wdl) | Stand-alone WDL workflows |   

---  

## Storage

All sample-level input data is stored in a secure Google Cloud bucket. Note that permissions must be granted for bucket access.  

The bucket is organized as follows:
* One directory per cohort
* Each cohort has the following subdirectories:
    * `gatk-hc/reblocked` : reblocked GATK-HC gVCFs and indexes
    * `gatk-sv` : evidence and metrics files collected by GATK-SV  
        * `gatk-sv/coverage` : coverage counts files  
        * `gatk-sv/metrics` : per-sample metrics generated during GATK-SV module 01  
        * `gatk-sv/pesr` : PE/SR metadata files  
    * `manta` : raw Manta VCFs and indexes  
    * `melt` : raw MELT VCFs and indexes  
  * `wham` : raw Wham VCFs and indexes  

_Note: all raw gVCFs previously hosted in `gatk-hc/` were deleted on Feb 7, 2024, but the nested directory structure of `gatk-hc/reblocked/` was preserved for legacy code compatability_  
