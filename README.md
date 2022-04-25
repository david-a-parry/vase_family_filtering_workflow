# segregation_workflow

[Snakemake](snakemake.github.io) workflow for annotating and filtering variant data in small families according to different inheritance patterns using [VASE](https://github.com/david-a-parry/vase).

## Introduction

[VASE](https://github.com/david-a-parry/vase) is a flexible tool for filtering variant data with a large number of filtering options, generally geared towards the analysis of variant data from small families with rare genetic disease. It has a potentially overwhelming number of options so this workflow aims to provide a general framework for the filtering of variant data which can be easily configuired using a YAML configuration file.

## Setup

### Snakemake

You will need to initialise a conda environment with Snakemake and pandas
available to run this workflow, ideally using [mamba](https://github.com/mamba-org/mamba).

To create the environment (only has to be run once):

    $ mamba create --name snakemake  'snakemake-minimal>=5.24.1' 'pandas>=1.1'

And then to activate this environment when you want to run this workflow:

    $ conda activate snakemake

The necessary environements for the different stages of this workflow  will then
be created when you run the workflow (see below).

## Running the workflow.

Before running the workflow you need to edit the [config/config.yaml](https://github.com/david-a-parry/vase_family_filtering_workflow/blob/main/config/config.yaml) file with your project specific settings. The required input files are:

* vcf: a [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) annotated [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) file containing variant calls for your samples of interest. This could potentially be generated using [this workflow](https://github.com/david-a-parry/dna-seq-gatk-variant-calling)
* ped: a [PED](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972) file detailing samples, their familial relationships and affected status

An example config file is provided in [config/config.yaml](https://github.com/david-a-parry/vase_family_filtering_workflow/blob/main/config/config.yaml) and a schema file detailing the available arguments can be found at [workflow/schemas/config.schema.yaml](https://github.com/david-a-parry/vase_family_filtering_workflow/blob/main/workflow/schemas/config.schema.yaml).

The workflow can then be run as you would any other snakemake workflow, for example:

    $ snakemake --use-conda

## Results

The `results/report` directory contains the final variant report in Excel (XLSX) and JSON format. The spreadsheet contains one worksheet per family with details of all variants meeting the filtering criteria. The JSON file has a similar structure for ease of use in any further programmatic analysis.

Several intermediate files are also provided for tweaking and rerunning analyses.

The `results/vase_annotated` directory contains a VCF annotated with information from any gnomAD, dbSNP, CADD and SpliceAI files used (`results/vase_annotated/vase_anno.vcf.gz`). The annotation steps are the most time-consuming due to the huge volume of data from these databases so this file is useful to have should you want to tweak any of your filtering parameters and rerun the remainder of the workflow. This directory also contains files with variants that could not be found in your local CADD/SpliceAI databases should you wish to manually score these variants for future analyses.

The `results/vase_filtered` directory contains VCF format files of all the variants passing the inheritance models specified in your config file.

## Issues

Because VASE is installed to its conda environment via pip, an issue can arise where conda creates a shebang line in the VASE executable that is too long, causing commands to fail. The shebang line points to the location of the python executable in the conda environment created by snakemake within the project folder. As such, the simplest fix that I am aware of is to run this workflow from a location with a relatively short path.


## Author

David A. Parry, University of Edinburgh



