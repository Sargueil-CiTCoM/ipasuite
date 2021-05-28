from snakemake.utils import validate
import pandas as pd
import os

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

unindexed_samples = pd.read_csv(config["samples"], sep="\t")#.set_index("sample", drop=False)
samples = unindexed_samples.set_index(["rna_id", "probe", "temperature", "magnesium", "replicate"],drop=True)
#samples = samples[samples.ref.notnull()]
#samples.index.names = ["sample_id"]

# TODO : Activate before prod
#validate(samples, schema="../schemas/samples.schema.yaml")

def get_raw_probe_input(wildcards):
    print(samples.index)
    file = samples.loc[(wildcards.rna_id, wildcards.probe, int(wildcards.temperature), wildcards.magnesium, int(wildcards.replicate))].dropna()
    print(file)
    return os.path.join(config["rawdata"]["path_prefix"], file["probe_file"]) 

def get_raw_control_input(wildcards):
    file = samples.loc[(wildcards.rna_id, wildcards.probe, int(wildcards.temperature), wildcards.magnesium, int(wildcards.replicate))].dropna()
    return os.path.join(config["rawdata"]["path_prefix"] + file["control_file"])

def get_qushape_refseq(wildcards):
    return None

def get_qushape_refproj(wildcards):
    return None

def get_all_raw_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        print(row)
        sample = ("resources/{fluodir}/{rna_id}" + config["format"]["condition"] + "_{replicate}.{fluoext}.tsv").format(fluodir = config["folders"][config["rawdata"]["type"]],
                fluoext=config["rawdata"]["type"],
                **row)
        control = ("resources/{fluodir}/{rna_id}" + config["format"]["control_condition"] + "_{replicate}.{fluoext}.tsv").format(control = config["rawdata"]["control"], fluodir =
                config["folders"][config["rawdata"]["type"]], fluoext=config["rawdata"]["type"],
                **row)
        
        outputs.append(sample)
        outputs.append(control)
    return outputs
        
#def get_ceq8000_input():
#    if config.rawdata["type"] == "fluo-ceq8000":
        

def get_final_outputs():
    exppath = "resources/{fluo-ceq8000}ceq8000/{folder}/{sample}.tsv"
    outputs = []
    print(samples)
    for row in samples.itertuples(name="Sample"):
        folder = row.folder
        sample = exppath.format(folder=folder,sample=os.path.splitext(row.sample_filename)[0])
        control = exppath.format(folder=folder,sample=os.path.splitext(row.control_filename)[0])
    return outputs
