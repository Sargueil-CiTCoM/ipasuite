from snakemake.utils import validate
import pandas as pd
import os

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

def get_indexes(replicates=True):
    indexes = ["rna_id"]
    indexes.extend(config["condition_names"])
    if replicates:
        indexes.append("replicate")
    return indexes


unindexed_samples = pd.read_csv(config["samples"], sep="\t")#.set_index("sample", drop=False)
samples = unindexed_samples.set_index(get_indexes(),drop=True)
samples_replicates = unindexed_samples.set_index(get_indexes(replicates=False),drop=True)

pool_ids = [pool["id"] for pool in config["ipanemap"]["pools"]]
#samples = samples[samples.ref.notnull()]
#samples.index.names = ["sample_id"]

# TODO : Activate before prod
#validate(samples, schema="../schemas/samples.schema.yaml")

def get_sample(wildcards, all_replicates=False):
    sval = [wildcards.rna_id]
    sval.extend([wildcards[cond] for cond in config["condition_names"]])

    if all_replicates:
        return samples_replicates.loc[tuple(sval)]

    sval.append(int(wildcards.replicate))
    return samples.loc[tuple(sval)]

def construct_path(step, control = False, results_dir = True, ext = None, replicate = True):
    cond = "_" + CONDITION if not control else expand("_" + CTRL_CONDITION, control=CONTROL, allow_missing = True)[0]
    basedir = "results" if results_dir else "resources"
    replicate = "_{replicate}" if replicate else ""
    extension = ".{step}.tsv" if ext is None else ext
    print(expand(basedir + "/{folder}/{rna_id}" + cond + replicate + extension, folder = FOLDERS[step], step=step, allow_missing=True))
    return expand(basedir + "/{folder}/{rna_id}" + cond + replicate + extension, folder = FOLDERS[step], step=step, allow_missing=True)


def get_raw_probe_input(wildcards):
    sample = get_sample(wildcards) 
    return os.path.join(config["rawdata"]["path_prefix"], sample["probe_file"]) 

def get_raw_control_input(wildcards):
    sample = get_sample(wildcards)
    return os.path.join(config["rawdata"]["path_prefix"] + sample["control_file"])
    

def get_qushape_refseq(wildcards):
    sample = get_sample(wildcards)
    path = config["sequences"][wildcards.rna_id]
    if os.path.exists(path):
        return path
    return []

def get_qushape_refproj(wildcards):
    sample = get_sample(wildcards)
    path = None
    if isinstance(sample["reference_qushape_file"], str):
        path = os.path.join(config["rawdata"]["path_prefix"] + sample["reference_qushape_file"])
    if not path is None and os.path.exists(path):
        return path
    return []

def get_replicates(wildcards, qushape_analysed = False):
    replicates = get_sample(wildcards, all_replicates=True)
    if qushape_analysed:
        return replicates.loc[replicates["qushape_analysed"] == "yes"]["replicate"]
    else:
        return replicates["replicate"]

def get_ipanemap_inputs(wildcards):
    inputs = []
    for pool in config["ipanemap"]["pools"]:
        if pool["id"] == wildcards.pool_id:
            for cond in pool["conditions"]:
                inputs.extend(expand(construct_path("aggreact-ipanemap", replicate=False, ext=".txt"),
                        rna_id=pool["rna_id"], **cond))
            break;
    return inputs

def get_all_raw_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        sample = ("resources/{fluodir}/{rna_id}" + config["format"]["condition"] + "_{replicate}.{fluoext}.tsv").format(fluodir = config["folders"][config["rawdata"]["type"]],
                fluoext=config["rawdata"]["type"],
                **row)
        control = ("resources/{fluodir}/{rna_id}" + config["format"]["control_condition"] + "_{replicate}.{fluoext}.tsv").format(control = config["rawdata"]["control"], fluodir =
                config["folders"][config["rawdata"]["type"]], fluoext=config["rawdata"]["type"],
                **row)
        
        outputs.append(sample)
        outputs.append(control)
    return outputs

def get_all_qushape_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        sample = ("results/{folder}/{rna_id}_" + config["format"]["condition"] + "_{replicate}.qushape").format(folder = config["folders"]["qushape"], **row)
        outputs.append(sample)
    return outputs

def get_all_reactivity_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        sample = ("results/{folder}/{rna_id}_" + config["format"]["condition"] + "_{replicate}.{step}.tsv").format(folder = config["folders"]["reactivity"],step="reactivity", **row)
        if row["qushape_analysed"] == 'yes':
            outputs.append(sample)
    return outputs
        
def get_all_aggreact_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        sample = expand(construct_path(step="aggreact", replicate=False), **row)
        sampleipan = expand(construct_path("aggreact-ipanemap", replicate=False, ext=".txt"),**row)
        if row["qushape_analysed"] == 'yes':
            outputs.append(sample)
            outputs.append(sampleipan)
        else:
            sample = expand(construct_path("qushape", ext=".qushape"), **row)
            outputs.append(sample)
    return outputs

## Not parallel, to be improved
def get_all_structure_outputs(wildcards):
    outputs = []
    for pool in config["ipanemap"]["pools"]:
        checkpoint_output = checkpoints.ipanemap.get(rna_id=pool["rna_id"], pool_id=pool["id"]).output[0]
        glob = glob_wildcards(os.path.join(checkpoint_output, "{rna_id}_pool_{pool_id}_optimal_{idx, \d+}.dbn"))
        outputs.extend(expand("results/{folder}/{rna_id}_pool_{pool_id}_{idx}.dbn",
                folder=config["folders"]["structure"], rna_id=pool["rna_id"], pool_id=pool["id"], idx=glob.idx))            

    return outputs 


#def get_final_outputs():
#    exppath = "resources/{fluo-ceq8000}ceq8000/{folder}/{sample}.tsv"
#    outputs = []
#    print(samples)
#    for row in samples.itertuples(name="Sample"):
#        folder = row.folder
#        sample = exppath.format(folder=folder,sample=os.path.splitext(row.sample_filename)[0])
#        control = exppath.format(folder=folder,sample=os.path.splitext(row.control_filename)[0])
#    return outputs
