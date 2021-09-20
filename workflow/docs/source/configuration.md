# Configuration 

The workflow contains 2 configuration files

- `config/config.yml` Which contains general configuration for the workflow
- `config/samples.tsv` That contains informations about samples and replicates

## Quickstart Workflow configuration 

`config.yaml` contains all informations necessary to execute the workflow. 
In this section we will discuss the main configurations


### 
sequences
: 


## Declare samples 

`.tsv` file is a tabular file format that you can open either with a text editor or an spreadsheet program

This file contains each experiment with related information. One line per experiment.

A `config/samples.tpl.tsv` is available, and can be used as a model for your project

### Mandatory columns

id (string)
: a unique number to identity experiment

rna_id (string)
: The identifier for RNA fragment used in this experiment, as declare in the `sequences` section of `config.yaml`

date (date)
: Date of the experiment. use ISO 8601 format YYYY-MM-DD 

replicate (integer)
: Replicate number for the given experiment


### Conditions columns
For each type of experimental condition, you must declare it in the `condition_names` of `config.yaml` file the name declared in the config file must be the same as the on in `samples.tsv`
In order to generate unambiguous file name, you must also add the conditions in the `format` section




#### Examples :

probe (string)
: Which probe was used in this sample (1M7, BzCN, NMIA, etc.)

temperature (int)
: At which temperature was made the probing

magnesium (string)
: Did the sample buffer contained Magnesium during probing (Mg / noMg)?

interaction (string)
: What other molecule/RNA was present with the probed RNA during the probing step.

### Optional columns

probe_file (relative file path)
: relative path to the external file fluo-ce file corresponding to this shape sample. if this field is filled in and no file is present for this sample in `resources/1.1-fluo-ce` snakemake will try to import this file. `rawdata->path_prefix` (from `config.yaml`) is prefixed to the content of the to construct fullpath.

control_file (relative file path)
: relative path to the external file fluo-ce file corresponding to the control file of this shape sample. If this field is filled in and no file is present for this sample in `resources/1.1-fluo-ce` snakemake will try to import this file. `rawdata->path_prefix` (from `config.yaml`) is prefixed to the content of the to construct fullpath.

reference_qushape_file (file path)
: Reference QuShape project : will be used in QuShape to pre-generate peak calling and alignment

### Other columns
You can add as many columns to your `samples.tsv` as you wish, to help you classify and caracterize your data. Each column must have un unique name.


#### Examples

experimenter
: Who performed the sample acquisition
