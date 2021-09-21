from snakemake.utils import validate
import pandas as pd
import os

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

def get_indexes(replicates=True):
    indexes = ["rna_id"]
    indexes.extend(config["conditions"])
    if replicates:
        indexes.append("replicate")
    return indexes

condition_types = { name : 'str' for name in config['conditions']}

unindexed_samples = pd.read_csv(config["samples"], sep="\t", dtype = condition_types)#.set_index("sample", drop=False)
samples = unindexed_samples.set_index(get_indexes(),drop=True)
samples_replicates = unindexed_samples.set_index(get_indexes(replicates=False),drop=True)

pool_ids = [pool["id"] for pool in config["ipanemap"]["pools"]]
#samples = samples[samples.ref.notnull()]
#samples.index.names = ["sample_id"]

# TODO : Activate before prod
#validate(samples, schema="../schemas/samples.schema.yaml")


def construct_list_param(config_category, param_name):
    if param_name in config_category and len(config_category[param_name]) > 0:
        return f"--{param_name}='{str(config_category[param_name])}'"
    return ""

def construct_param(config_category, param_name):
    if param_name in config_category:
        return f"--{param_name}='{str(config_category[param_name])}'"
    return ""

## Relative to tools/normalize_reactivity
def construct_normcol():
    arg = ""
    if "norm_column" in config["aggregate"]:
        arg = f"--normcol={config['aggregate']['normcol']}"
    elif "norm_method" in config["aggregate"]:
        if "simple":
            arg =  "--normcol=simple_norm_reactivity"
        elif "interquartile":
            arg = "--normcol=interquartile_norm_reactivity"
    return arg


def get_sample(wildcards, all_replicates=False):
    sval = [wildcards.rna_id]
    sval.extend([wildcards[cond] for cond in config["conditions"]])

    if all_replicates:
        return samples_replicates.loc[tuple(sval)]

    sval.append(int(wildcards.replicate))
    return samples.loc[tuple(sval)]

def construct_path(step, control = False, results_dir = True, ext = None, replicate = True, log_dir = False):
    cond = f'_{CONDITION}' if not control else expand(f'_{CTRL_CONDITION}', control= CONTROL, allow_missing = True)[0]
    basedir = "logs" if log_dir else ("results" if results_dir else "resources")
    replicate = "_{replicate}" if replicate else ""
    extension = ".{step}.tsv" if ext is None else ext
    pid = '{folder}/{rna_id}' if not log_dir else '{step}-{rna_id}'
    path = expand(f'{basedir}/{pid}{cond}{replicate}{extension}', folder = FOLDERS[step], step=step, allow_missing=True)
    #print(path)
    return path

def get_external_qushape(wildcards):
    sample = get_sample(wildcards)
    return os.path.join(config["rawdata"]["path_prefix"], sample["qushape_file"])

def get_raw_probe_input(wildcards):
    sample = get_sample(wildcards)
    return os.path.join(config["rawdata"]["path_prefix"], sample["probe_file"])

def get_raw_control_input(wildcards):
    sample = get_sample(wildcards)
    return os.path.join(config["rawdata"]["path_prefix"], sample["control_file"])


def get_refseq(wildcards, all_replicates = False):
    sample = get_sample(wildcards, all_replicates)
    path = config["sequences"][wildcards.rna_id]
    if os.path.exists(path):
        return path
    return []

def get_qushape_refproj(wildcards):
    sample = get_sample(wildcards)
    path = None
    if isinstance(sample["reference_qushape_file"], str):
        path = os.path.join(config["rawdata"]["path_prefix"], sample["reference_qushape_file"])
    if not path is None and os.path.exists(path):
        print(path)
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
        sample = construct_path(results_dir=False, step=RAW_DATA_TYPE)[0].format(**row)
        control = construct_path(results_dir=False, step=RAW_DATA_TYPE, control=True)[0].format(**row)
        outputs.append(sample)
        outputs.append(control)
    return outputs

def get_all_qushape_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        sample = construct_path(step="qushape", ext=".qushape")[0].format(**row)
        outputs.append(sample)
    return outputs

def get_all_reactivity_outputs():
    outputs = []
    for idx,row in samples.reset_index().iterrows():
        sample = ("results/{folder}/{rna_id}_" + 
                  config["format"]["condition"] + 
                  "_{replicate}.{step}.tsv").format(
                          folder = config["folders"]["reactivity"],
                          step="reactivity",
                          **row)
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

def get_all_vienna_outputs(wildcards):
    outputs = []
    for pool in config["ipanemap"]["pools"]:
        checkpoint_output = checkpoints.ipanemap.get(rna_id=pool["rna_id"], pool_id=pool["id"]).output[0]
        glob = glob_wildcards(os.path.join(checkpoint_output, "{rna_id}_pool_{pool_id}_optimal_{idx, \d+}.dbn"))
        outputs.extend(expand("results/{folder}/{rna_id}_pool_{pool_id}_{idx}.varna",
            folder=config["folders"]["varna"], rna_id=pool["rna_id"], pool_id=pool["id"], idx=glob.idx))
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
