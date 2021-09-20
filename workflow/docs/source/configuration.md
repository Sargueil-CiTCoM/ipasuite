# Configuration 

The workflow contains 2 configuration files

- `config/config.yml` Which contains general configuration for the workflow
- `config/samples.tsv` That contains informations about samples and replicates

## Configure Pipeline

`config.yaml` contains all informations necessary to execute the workflow

[TODO]

## Declare samples 

`.tsv` file is a tabular file format that you can open either with a text editor or an spreadsheet program

This file contains each experiment with related information. One line per experiment.

A `config/samples.tpl.tsv` is available, and can be used as a model for your project


### Mandatory columns

id (string)
: a unique number to identity experiment

rna_id (string)
: The identifier for RNA fragment used in this experiment

probe (string)
: Which probe has been use for this experiment (1M7, NMIA, BzCN, etc.)

date (date)
: Date of the experiment. use ISO 8601 format YYYY-MM-DD 

replicate (integer)
: Replicate number for the given experiment


### Special columns

condition column (string)
: One column per experimental condition ( temperature, magnesium, probing delay, folding buffer) The value in those column will be use to construct filename 

### Optional columns

probe_file (relative file path)
: relative path to the external file fluo-ce file corresponding to this shape sample. if this field is filled in and no file is present for this sample in `resources/1.1-fluo-ce` snakemake will try to import this file. `rawdata->path_prefix` (from `config.yaml`) is prefixed to the content of the to construct fullpath.

control_file (relative file path)
: relative path to the external file fluo-ce file corresponding to the control file of this shape sample. If this field is filled in and no file is present for this sample in `resources/1.1-fluo-ce` snakemake will try to import this file. `rawdata->path_prefix` (from `config.yaml`) is prefixed to the content of the to construct fullpath.

reference_qushape_file (file path)
: Reference QuShape project : will be used in QuShape to pre-generate peak calling and alignment

### Other columns
You can add as many columns to your `samples.tsv` as you wish  

All other column are optionaluseful for the experimenter, but not mandatory to perform data treatement.

### Example columns :

| id  | rna_id   | experimenter | probe | temperature | magnesium | folding_buffer | probing_delay | primer | rna_begin | rna_end | date       | replicate | qushape_analysed | reference_qushape_file                  | reference_sequence_file | probe_file                              | control_file                                    |
| --- | -------- | ------------ | ----- | ----------- | --------- | -------------- | ------------- | ------ | --------- | ------- | ---------- | --------- | ---------------- | --------------------------------------- | ----------------------- | --------------------------------------- | ----------------------------------------------- |
| 1   | didymium |              | 1M7   | 37C         | Mg        |                |               |        |           |         | 2019-08-01 | 1         | yes              | 2-qushape/didymium_1M7_37C_Mg_3.qushape | didymium.fa             | input/didymium_1M7_37C_Mg_1.fluo-ce.tsv | input/didymium_DMSO_of_1M7_37C_Mg_1.fluo-ce.tsv |









