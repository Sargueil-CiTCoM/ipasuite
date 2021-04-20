from snakemake.utils import validate
import pandas as pd
import os

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t")#.set_index("sample", drop=False)
samples = samples[samples.ref.notnull()]
#samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

def get_final_outputs():
    exppath = "resources/ceq8000/{folder}/{sample_filename}.tsv"
    outputs = expand(exppath, sample_filename=samples["sample_filename"].apply(lambda x: os.path.splitext(x)[0]), folder=samples["folder"])
    outputs.append(expand(exppath, sample_filename=samples["control_filename"].apply(lambda x: os.path.splitext(x)[0]), folder=samples["folder"]))

    return outputs
