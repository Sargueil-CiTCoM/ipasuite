from snakemake.utils import validate
import pandas as pd
import os


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### load config and sample sheets #####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")


CONDITION = config["format"]["condition"]
CTRL_CONDITION = config["format"]["control_condition"]
RAW_DATA_TYPE = config["rawdata"]["type"]
CONTROL = config["rawdata"]["control"]
FOLDERS = config["folders"]
RESULTS_DIR = config["results_dir"] if "results_dir" in config else "results"
RESOURCES_DIR = config["resources_dir"] if "resources_dir" in config else "resources"
MESSAGE = config["format"]["message"]
CNORM = config["normalization"]

def get_indexes(replicates_in_index=True):
    indexes = ["rna_id"]
    indexes.extend(config["conditions"])
    if replicates_in_index:
        indexes.append("replicate")
    return indexes


condition_types = {name: "str" for name in config["conditions"]}
condition_types['id'] = "str"
condition_types['rna_id'] = "str"
condition_types['replicate'] = "str"

unindexed_samples = pd.read_csv(
    config["samples"], sep="\t", dtype=condition_types
)  # .set_index("sample", drop=False)

unindexed_samples = unindexed_samples.loc[ \
    (unindexed_samples["discard"] != "yes") & \
    (unindexed_samples["discard"] != True)]

validate(unindexed_samples, schema="../schemas/samples.schema.yaml",
        set_default=True)

if config['qushape']['use_subsequence']:
    unindexed_samples['rt_begin_pos'] = unindexed_samples['rt_begin_pos'].astype(int)
    unindexed_samples['rt_end_pos'] = unindexed_samples['rt_end_pos'].astype(int)

samples = unindexed_samples.set_index(get_indexes(), drop=True)
index_count = pd.Index(samples.index).value_counts()
if len(index_count[index_count > 1]) > 0:
    print("Duplicated entry in samples file")
    print(f"Entries: {index_count[index_count > 1]}")
    raise Exception("Duplicated entry")

samples_no_rep_index = unindexed_samples.set_index(
    get_indexes(replicates_in_index=False), drop=True
)

pool_ids = [pool["id"] for pool in config["ipanemap"]["pools"]]
# samples = samples[samples.ref.notnull()]
# samples.index.names = ["sample_id"]



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
            arg = "--normcol=simple_norm_reactivity"
        elif "interquartile":
            arg = "--normcol=interquartile_norm_reactivity"
    return arg


def get_sample(wildcards, list_replicates=False):
    sval = [wildcards.rna_id]
    sval.extend([wildcards[cond] for cond in config["conditions"]])

    if list_replicates:
        return samples_no_rep_index.loc[[tuple(sval)]]
    
    sval.append(wildcards.replicate)
    return samples.loc[[tuple(sval)]]


def construct_path(
    step, control=False, results_dir=True, ext=None, show_replicate=True,
    log_dir=False, figure=False, split_seq=False):
    cond = (
        f"_{CONDITION}"
        if not control
        else expand(f"_{CTRL_CONDITION}", control=CONTROL, allow_missing=True)[0]
    )
    figdir = "figures/" if figure else ""
    basedir = "logs" if log_dir else (RESULTS_DIR if results_dir else RESOURCES_DIR)
    replicate = "_{replicate}" if show_replicate else ""
    extension = ".{step}.tsv" if ext is None else ext
    pid = "{folder}/{rna_id}" if not log_dir else "{step}-{rna_id}"
    ssplit = "_seq{rt_end_pos}-{rt_begin_pos}" if split_seq \
             and config['qushape']['use_subsequence'] \
             else "" 

    path = expand(
        f"{basedir}/{figdir}{pid}{cond}{replicate}{ssplit}{extension}",
        folder=FOLDERS[step],
        step=step,
        allow_missing=True,
    )
    # print(path)
    return path

def aggregate_input_type():
    if config["qushape"]["use_subsequence"]:
        return "alignnormreact"
    else:
        return "normreact"


def get_external_qushape(wildcards):
    sample = get_sample(wildcards).iloc[0]
    return os.path.join(config["rawdata"]["path_prefix"], sample["qushape_file"])


def get_raw_probe_input(wildcards):
    sample = get_sample(wildcards).iloc[0]
    return os.path.join(config["rawdata"]["path_prefix"], sample["probe_file"])


def get_raw_control_input(wildcards):
    sample = get_sample(wildcards).iloc[0]
    return os.path.join(config["rawdata"]["path_prefix"], sample["control_file"])

def get_align_begin(wildcards):
    sample = get_sample(wildcards).iloc[0]
    return sample["rt_end_pos"]

def get_align_end(wildcards):
    sample = get_sample(wildcards).iloc[0]
    return sample["rt_begin_pos"]

def get_align_reactivity_inputs(wildcards):
    sample = get_sample(wildcards).iloc[0]
    fa = f"{RESULTS_DIR}/{config['folders']['subseq']}/" \
         f"{wildcards.rna_id}_{sample['rt_end_pos']}-{sample['rt_begin_pos']}.fasta"
    norm = expand(construct_path("normreact", split_seq=True),
            rt_end_pos=sample['rt_end_pos'],rt_begin_pos=sample['rt_begin_pos'],
            **wildcards )
    return {"refseq": fa, "norm": norm} 


def get_subseq(wildcards, split_seq=False):
    sample = get_sample(wildcards, list_replicates=split_seq)
    fasta = config["sequences"][wildcards.rna_id]

    if split_seq and config["qushape"]["use_subsequence"]:
        return  f"{RESULTS_DIR}/{config['folders']['subseq']}/" \
                f"{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.fasta"

    if os.path.exists(fasta):
        return fasta
    return []

def get_refseq(wildcards):
    fasta = config["sequences"][wildcards.rna_id]
    if os.path.exists(fasta):
        return fasta
    return []



def get_qushape_refproj(wildcards):
    sample = get_sample(wildcards).iloc[0]
    path = None
    if isinstance(sample["reference_qushape_file"], str):
        path = os.path.join(
            config["rawdata"]["path_prefix"], sample["reference_qushape_file"]
        )

    if path is not None and os.path.exists(path):
        return path
    return []


def get_replicate_list(wildcards):
    replicates = get_sample(wildcards, list_replicates=True)
    return replicates["replicate"]


def get_ipanemap_pool_inputs(pool):
    inputs = []
    for cond in pool["conditions"]:
        inputs.extend(
            expand(
                construct_path(
                    "aggreact-ipanemap", show_replicate=False, ext=".txt"
                ),
                rna_id=pool["rna_id"],
                **cond,
            )
        )
    return inputs


def get_ipanemap_inputs(wildcards):
    for pool in config["ipanemap"]["pools"]:
        if pool["id"] == wildcards.pool_id:
            return get_ipanemap_pool_inputs(pool) 
    return []


def get_all_raw_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(results_dir=False, step=RAW_DATA_TYPE)[0].format(**row)
        control = construct_path(results_dir=False, step=RAW_DATA_TYPE, control=True)[
            0
        ].format(**row)
        outputs.append(sample)
        outputs.append(control)
    return outputs


def get_all_qushape_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(step="qushape", ext=".qushape", split_seq=True)[0].format(**row)
        outputs.append(sample)
    return outputs


def get_all_reactivity_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = (
            f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_{config['format']['condition']}_{{replicate}}.{{step}}.tsv"
        ).format(folder=config["folders"]["reactivity"], step="reactivity", **row)
        outputs.append(sample)
    return outputs


def get_all_aggreact_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = expand(construct_path(step="aggreact", show_replicate=False), **row)
        sampleipan = expand(
            construct_path("aggreact-ipanemap", show_replicate=False, ext=".txt"), **row
        )
        outputs.append(sample)
        outputs.append(sampleipan)
#        else:
#            sample = expand(construct_path("qushape", ext=".qushape"), **row)
#            outputs.append(sample)
    return outputs


## Not parallel, to be improved
def get_all_structure_outputs(wildcards):
    outputs = []
    for pool in config["ipanemap"]["pools"]:
        checkpoint_output = checkpoints.ipanemap.get(
            rna_id=pool["rna_id"], pool_id=pool["id"]
        ).output[0]
        glob = glob_wildcards(
            os.path.join(
                checkpoint_output, "{rna_id}_pool_{pool_id}_optimal_{idx, \d+}.dbn"
            )
        )
        outputs.extend(
            expand(
                f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx}}.dbn",
                folder=config["folders"]["structure"],
                rna_id=pool["rna_id"],
                pool_id=pool["id"],
                idx=glob.idx,
            )
        )

    return outputs

def get_varna_pool_concat_inputs(wildcards):
    inputs = []
    for pool in config["ipanemap"]["pools"]:
        if pool["id"] == wildcards.pool_id:
            for cond in pool["conditions"]:
                inputs.extend(
                    expand(
                        f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx}}_cond_{CONDITION}.varna",
                        folder=config["folders"]["varna"],
                        **wildcards,
                        **cond
                    )
                )
            return inputs
    return inputs


def get_all_varna_by_condition_outputs(wildcards):
    outputs = []
    for pool in config["ipanemap"]["pools"]:
        checkpoint_output = checkpoints.ipanemap.get(
            rna_id=pool["rna_id"], pool_id=pool["id"]
        ).output[0]
        glob = glob_wildcards(
            os.path.join(
                checkpoint_output, "{rna_id}_pool_{pool_id}_optimal_{idx, \d+}.dbn"
            )
        )

        for cond in pool["conditions"]:
            outputs.extend(
                expand(
                    f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx}}_cond_{CONDITION}.varna",
                    folder=config["folders"]["varna"],
                    rna_id=pool["rna_id"],
                    pool_id=pool["id"],
                    idx=glob.idx,
                    **cond
                )
            )
#        outputs.extend(
#            expand(
#                "{RESULTS_DIR}/{folder}/{rna_id}_pool_{pool_id}_{idx}.svg",
#                folder=config["folders"]["varna"],
#                rna_id=pool["rna_id"],
#                pool_id=pool["id"],
#                idx=glob.idx,
#            )
#        )
    return outputs

def get_all_varna_pool_concat_outputs(wildcards):
    outputs = []
    for pool in config["ipanemap"]["pools"]:
        checkpoint_output = checkpoints.ipanemap.get(
            rna_id=pool["rna_id"], pool_id=pool["id"]
        ).output[0]
        glob = glob_wildcards(
            os.path.join(
                checkpoint_output, "{rna_id}_pool_{pool_id}_optimal_{idx, \d+}.dbn"
            )
        )

        outputs.extend(
            expand(
                f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_pool_{{pool_id}}_{{idx}}.varna",
                folder=config["folders"]["varna"],
                rna_id=pool["rna_id"],
                pool_id=pool["id"],
                idx=glob.idx,
            )
        )
#        outputs.extend(
#            expand(
#                "{RESULTS_DIR}/{folder}/{rna_id}_pool_{pool_id}_{idx}.svg",
#                folder=config["folders"]["varna"],
#                rna_id=pool["rna_id"],
#                pool_id=pool["id"],
#                idx=glob.idx,
#            )
#        )
    return outputs

# def get_final_outputs():
#    exppath = "resources/{fluo-ceq8000}ceq8000/{folder}/{sample}.tsv"
#    outputs = []
#    print(samples)
#    for row in samples.itertuples(name="Sample"):
#        folder = row.folder
#        sample = exppath.format(folder=folder,sample=os.path.splitext(row.sample_filename)[0])
#        control = exppath.format(folder=folder,sample=os.path.splitext(row.control_filename)[0])
#    return outputs
