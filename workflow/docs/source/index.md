# Welcome to Shape-CE SnakeMake Pipeline documentation!

This documentation will walk you through the step needed to execute the shape-ce pipeline

## Presentation

Shape-CE workflow intend to provide some automation in the data treatment of SHAPE Capillary Electrophorese.

It rely on : 

- QuShape
- IPANEMAP
- RNAFold
- VARNA. 
- Custom scripts for file conversion, reactivity normalization and aggregation.

The workflow will enable you to generate structure data for a RNA fragment analysed using SHAPE with a arbitrary set of conditions (Temperature, Magnesium, Probes.)
Data will be organized using a `config/samples.tsv` and configured using `config/config.yaml` , and be feeded to IPANEMAP with a fully configured set of conditions.

### Use this pipeline for your project

- Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
- Clone the newly created repository to your local system, into the place where you want to perform the data analysis.

# Preparation

## Operating System

The pipeline should work on any system suporting Conda. 
However, some packages (ViennaRNA, scikit-bio) are not available in Windows using conda, you will have install them manually.

Tested Operating system :
- Archlinux
- Debian 11 (Bulleyes)
- Ubuntu 16.4 to 20.4
- Windows 10



## Install Conda environnement

- Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

- Open a terminal configured with conda and go to the root of the project

- Prepare Snakemake conda environnement
  
  ```
  conda env create --name snakemake --file workflow/envs/snakemake.yml
  ```

## Prepare workflow

### File organization

```
├── config # Contains configurations files
│   ├── config.yaml
│   └── samples.tsv
├── resources # Contains input files
│   ├── 1.1-fluo-ceq8000 # SHAPE data from ceq8015
│   ├── 1.2-fluo-ce # SHAPE data in a tsv format readable by QuShape
│   ├── sequence1.fa # Reference sequences of the studied RNA
│   ├── ... 
│   └── sequenceN.fa
├── results # Intermediary and analysied data
│   ├── 2-qushape
│   ├── 3.1-reactivity
│   ├── 3.2-normreact
│   ├── 4.1-aggreact
│   ├── 4.2-aggreact-ipanemap
│   ├── 5.1-ipanemap-config
│   ├── 5.2-ipanemap-out
│   ├── 5.3-structure
│   └── 5.4-varna 
└── workflow # Scripts used for data treatement
```

The workflow contains 2 configuration files

- `config/config.yml` Which contains general configuration for the workflow
- `config/samples.tsv` That contains informations about samples and replicates

### Sample configuration file

`.tsv` file is a tabular file format that you can open either with a text editor or an spreadsheet program

This file contains each experiment with related information. One line per experiment.

Mandatory columns  :

| name                    | description                                                                                                                                                                                                         |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| id                      | a unique number to identity experiment                                                                                                                                                                              |
| rna_id                  | The identifier for RNA fragment used in this experiment                                                                                                                                                             |
| probe                   | Which probe has been use for this experiment (1M7, NMIA, BzCN, etc.)                                                                                                                                                |
| [conditions]            | One column per experimental condition ( temperature, magnesium, probing delay, folding buffer), column values must be string. Do not use boolean value. The value in those column will be use to construct filename |
| date                    | Date of the experiment. use ISO 8601 format YYYY-MM-DD                                                                                                                                                              |
| replicate               | Replicate number for the given experiment                                                                                                                                                                           |
| qushape_analysed        | (yes/no) indicate whether the qushape analysis as be perform or not. If **no** the pipeline will only perform step until qushape project generations. If **yes** analysis will be performed entierly                |
| probe_file              | path to the fluo-ce file corresponding to this shape experiment (1M7, NMIA, BzCN, etc.)                                                                                                                             |
| control_file            | path to the control fluo-ce file corresponding to this shape experiment (DMSO)                                                                                                                                      |
| reference_qushape_file  | Reference QuShape project : will be used in QuShape to pre-generate peak calling and alignment                                                                                                                      |
| reference_sequence_file | Reference sequence feeded to QuShape. if not filled in the global sequence from config.yaml will be used                                                                                                            |

 All other column are useful for the experimenter, but not mandatory to perform data treatement.

#### Example columns :

| id  | rna_id   | experimenter | probe | temperature | magnesium | folding_buffer | probing_delay | primer | rna_begin | rna_end | date       | replicate | qushape_analysed | reference_qushape_file                  | reference_sequence_file | probe_file                              | control_file                                    |
| --- | -------- | ------------ | ----- | ----------- | --------- | -------------- | ------------- | ------ | --------- | ------- | ---------- | --------- | ---------------- | --------------------------------------- | ----------------------- | --------------------------------------- | ----------------------------------------------- |
| 1   | didymium |              | 1M7   | 37C         | Mg        |                |               |        |           |         | 2019-08-01 | 1         | yes              | 2-qushape/didymium_1M7_37C_Mg_3.qushape | didymium.fa             | input/didymium_1M7_37C_Mg_1.fluo-ce.tsv | input/didymium_DMSO_of_1M7_37C_Mg_1.fluo-ce.tsv |

### Config.yaml

`config.yaml` contains all informations necessary to execute the workflow


``` yaml


```







```{toctree}
:maxdepth: 2
:caption: Contents
```

Indices and tables
==================

* {ref}`genindex`
Ê* {ref}`modindex`
* {ref}`search`
