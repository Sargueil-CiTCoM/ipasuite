# Run the pipeline

## General information

Once you configured your pipeline it is time to launch it.

Snakemake will try to go forward with all step available :

-   If input file are missing it will try to import them using `rawdata:prefix_path` from `config.yaml` concatenated with `probe_file` and `control_file` column from `samples.tsv`
-   If QuShape projects does not exist, it will created using sequencer data.
-   If QuShape projects exists, it will try extract reactivity and go on until structures are generated

Simplest way to run the pipeline is to launch `./shape-ce.sh` from the root of your project

```bash
conda activate snakemake 
# or
source activate snakemake

./shape-ce.sh
# or
snakemake -j8 --use-conda --keep-going
```

It will run `snakemake` with 8 threads, using conda and do as much treatment it can.
At the end of snakemake run, it will output logs if they contains informations.

In order to control when the pipeline will stop during the treatement, you can make it stop at intermediate steps :

To stop after generating QuShape project : `./shape-ce.sh all_qushape`
To stop after sample aggregation (before IPANEMAP) : `./shape-ce.sh all_aggregate`

## In practice 

1. Open an terminal / console
2. Got to your git repository
```bash
cd path/to/project/
```
3. run qushape projects generator 
```
./shape-ce.sh all_qushape
```
4. Open each new qushape project with QuShape and perform data treatment until reactivity step (no reactivity export needed)
5. run again 
```./shape-ce.sh```

You will find generated structure in `results/5.3-structure` and varna file in `results/5.4-varna` 
