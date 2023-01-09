from snakemake.utils import validate
import snakemake
import os
import json
import load_samples
from functools import partial

# #### load config and sample sheets #####

if "config" not in vars():
    from snakemake.utils import expand
    from snakemake.checkpoints import checkpoints
    from snakemake.io import glob_wildcards

    config = {}


validate(config, schema="../schemas/config.schema.yaml")


if "previous" in config["format"]:
    PREV_CONDITION = config["format"]["previous"]["condition"]
    PREV_CTRL_CONDITION = config["format"]["previous"]["control_condition"]
CONDITION = config["format"]["condition"]
CTRL_CONDITION = config["format"]["control_condition"]
RAW_DATA_TYPE = config["rawdata"]["type"]
CONTROL = config["rawdata"]["control"]
FOLDERS = config["folders"]
RESULTS_DIR = config["results_dir"] if "results_dir" in config else "results"
RESOURCES_DIR = config["resource_dir"] if "resource_dir" in config else "resources"
MESSAGE = config["format"]["message"]
CNORM = config["normalization"]

indexes = load_samples.get_indexes(config, replicates_in_index=True)
indexes_no_rep = load_samples.get_indexes(config, replicates_in_index=False)
unindexed_samples = load_samples.get_unindexed_samples(config)
samples = load_samples.get_samples(config)
samples_no_rep_index = load_samples.get_samples(config, replicate_in_index=False)

get_sample = partial(load_samples.get_sample, config)


pool_ids = [pool["id"] for pool in config["ipanemap"]["pools"]]
# samples = samples[samples.ref.notnull()]
# samples.index.names = ["sample_id"]


def get_reactive_nucleotides(wildcards):
    if wildcards.probe in config["probe"]:
        return construct_list_param(
            config["probe"][wildcards.probe], "reactive_nucleotides"
        )
    else:
        return construct_list_param(CNORM, "reactive_nucleotides")


def construct_dict_param(config_category, param_name):
    if param_name in config_category and len(config_category[param_name]) > 0:
        dictparam = json.dumps(dict(config_category[param_name]))
        return f"--{param_name}='{dictparam}'"
    return ""


def construct_list_param(config_category, param_name):
    if param_name in config_category and len(config_category[param_name]) > 0:
        return f"--{param_name}='{str(config_category[param_name])}'"
    return ""


def construct_param(config_category, param_name):
    if param_name in config_category:
        return f"--{param_name}='{str(config_category[param_name])}'"
    return ""


# Relative to tools/normalize_reactivity
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


def construct_path(
    step,
    control=False,
    results_dir=True,
    ext=None,
    show_replicate=True,
    log_dir=False,
    figure=False,
    split_seq=False,
    force_split_seq=False,
    previous=False,
    merged_conditions=False,
    source=None,
):
    if merged_conditions:
        cond = "_{conditions}"
    elif previous:
        cond = (
            f"_{PREV_CONDITION}"
            if not control
            else expand(f"_{PREV_CTRL_CONDITION}", control=CONTROL, allow_missing=True)[
                0
            ]
        )
    else:
        cond = (
            f"_{CONDITION}"
            if not control
            else expand(f"_{CTRL_CONDITION}", control=CONTROL, allow_missing=True)[0]
        )
    figdir = "figures/" if figure else ""
    basedir = (
        f"{RESULTS_DIR}/logs"
        if log_dir
        else (RESULTS_DIR if results_dir else RESOURCES_DIR)
    )
    replicate = "_{replicate}" if show_replicate else ""
    extension = ".{step}.tsv" if ext is None else ext
    pid = "{folder}/{rna_id}" if not log_dir else "{step}-{rna_id}"
    ssplit = (
        "_seq{rt_end_pos}-{rt_begin_pos}"
        if (split_seq and config["qushape"]["use_subsequence"]) or force_split_seq
        else ""
    )

    path = expand(
        f"{basedir}/{figdir}{pid}{cond}{replicate}{ssplit}{extension}",
        folder=FOLDERS[step],
        step=step,
        allow_missing=True,
    )
    # print(step)
    # print(source)
    # print(merged_conditions)
    # print(path)
    return path


def construct_full_condition_path(
    wildcards, step, results_dir=True, ext=None, previous=False
):
    path = construct_path(
        step=step,
        results_dir=results_dir,
        merged_conditions=False,
        previous=previous,
        ext=ext,
    )[0]
    # print(path)
    if previous:
        return path  # f"--full-previous-condition-format {path}"
    else:
        return path  # f"--full-condition-format {path}"


def aggregate_input_type():
    if config["qushape"]["use_subsequence"]:
        return "alignnormreact"
    else:
        return "normreact"


def get_ddntp_qushape(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]
    return sample["ddNTP"]


def get_external_qushape(wildcards):
    sample = get_sample(samples, wildcards)  # .iloc[0]
    return os.path.join(config["rawdata"]["path_prefix"], sample["qushape_file"])


def get_raw_probe_input(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]
    return os.path.join(config["rawdata"]["path_prefix"], sample["probe_file"])


def get_raw_control_input(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]
    return os.path.join(config["rawdata"]["path_prefix"], sample["control_file"])


def get_align_begin(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]
    return sample["rt_end_pos"]


def get_align_end(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]
    return sample["rt_begin_pos"]


def get_align_reactivity_inputs(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]

    fa = f"{config['sequences'][wildcards.rna_id]}"
    # fa = f"{RESULTS_DIR}/{config['folders']['subseq']}/" \
    # f"{wildcards.rna_id}_{sample['rt_end_pos']}-{sample['rt_begin_pos']}.fasta"
    norm = expand(
        construct_path("normreact", split_seq=True),
        rt_end_pos=sample["rt_end_pos"],
        rt_begin_pos=sample["rt_begin_pos"],
        **wildcards,
    )
    return {"refseq": fa, "norm": norm}


def get_subseq(wildcards, split_seq=False):
    # sample = get_sample(wildcards, list_replicates=split_seq)
    fasta = config["sequences"][wildcards.rna_id]

    if split_seq and config["qushape"]["use_subsequence"]:
        return (
            f"{RESULTS_DIR}/{config['folders']['subseq']}/"
            f"{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.fasta"
        )

    if os.path.exists(fasta):
        return fasta
    else:
        print(f"Missing FASTA {fasta}")
        raise
    return []


def get_refseq(wildcards):
    fasta = config["sequences"][wildcards.rna_id]
    if os.path.exists(fasta):
        return fasta
    return []


def get_qushape_refproj(wildcards):
    sample = get_sample(samples, wildcards).iloc[0]
    path = None
    if isinstance(sample["reference_qushape_file"], str):
        path = os.path.join(
            config["rawdata"]["path_prefix"], sample["reference_qushape_file"]
        )

    if path is not None and os.path.exists(path):
        return path
    return []


def get_replicate_list(wildcards):
    replicates = get_sample(samples_no_rep_index, wildcards)
    return replicates["replicate"]


def get_ipanemap_pool_inputs(pool):
    inputs = []
    for cond in pool["conditions"]:
        inputs.extend(
            expand(
                construct_path("aggreact-ipanemap", show_replicate=False, ext=".shape"),
                rna_id=pool["rna_id"],
                **cond,
            )
        )
    return inputs


def get_footprint_comp_input(comp, cond_name):
    return expand(
        construct_path("aggreact", show_replicate=False),
        rna_id=comp["rna_id"],
        **comp[cond_name],
    )


def get_footprint_inputs(wildcards):
    inputs = []
    for comp in config["footprint"]["compares"]:
        if comp["id"] == wildcards.foot_id:
            inputs.extend(get_footprint_comp_input(comp, "condition1"))
            inputs.extend(get_footprint_comp_input(comp, "condition2"))

    return inputs


def get_footprint_condition_name(wildcards, idx):
    for comp in config["footprint"]["compares"]:
        if comp["id"] == wildcards.foot_id:
            return config["format"]["condition"].format(**comp[f"condition{idx}"])

    return None


def get_ipanemap_inputs(wildcards):
    for pool in config["ipanemap"]["pools"]:
        if pool["id"] == wildcards.pool_id:
            return get_ipanemap_pool_inputs(pool)
    return []


def get_all_raw_rename_outputs(raw_data_type=RAW_DATA_TYPE):
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(results_dir=False, step=raw_data_type)[0].format(**row)

        control = construct_path(results_dir=False, step=raw_data_type, control=True,)[
            0
        ].format(**row)
        sample_src = construct_path(
            results_dir=False, step=raw_data_type, previous=True
        )[0].format(**row)
        control_src = construct_path(
            results_dir=False, step=raw_data_type, control=True, previous=True
        )[0].format(**row)
        if os.path.exists(sample_src):
            outputs.append(sample)
        if os.path.exists(control_src):
            outputs.append(control)
    # print(outputs)
    return outputs


def get_all_raw_addpositions_outputs(raw_data_type=RAW_DATA_TYPE):
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(
            results_dir=False, step=raw_data_type, force_split_seq=True
        )[0].format(**row)

        control = construct_path(
            results_dir=False, step=raw_data_type, control=True, force_split_seq=True
        )[0].format(**row)
        sample_src = construct_path(
            results_dir=False,
            step=raw_data_type,
            force_split_seq=False,
            split_seq=False,
        )[0].format(**row)
        control_src = construct_path(
            results_dir=False,
            step=raw_data_type,
            control=True,
            force_split_seq=False,
            split_seq=False,
        )[0].format(**row)
        if os.path.exists(sample_src):
            outputs.append(sample)
        if os.path.exists(control_src):
            outputs.append(control)
    return outputs


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


def get_all_qushape_rename_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(step="qushape", ext=".qushape",)[
            0
        ].format(**row)
        sample_src = construct_path(step="qushape", ext=".qushape", previous=True)[
            0
        ].format(**row)
        if os.path.exists(sample_src):
            outputs.append(sample)
    return outputs


def get_all_qushape_addpositions_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(
            step="qushape", ext=".qushape", split_seq=True, force_split_seq=True
        )[0].format(**row)
        sample_src = construct_path(
            step="qushape", ext=".qushape", split_seq=False, force_split_seq=False
        )[0].format(**row)
        if os.path.exists(sample_src):
            outputs.append(sample)
    return outputs

def get_all_qushapey_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(step="qushape", ext=".qushapey", split_seq=True)[
            0
        ].format(**row)
        outputs.append(sample)
    return outputs

def get_all_qushape_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(step="qushape", ext=".qushape", split_seq=True)[
            0
        ].format(**row)
        outputs.append(sample)
    return outputs


def get_all_reactivity_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = construct_path(step="reactivity", split_seq=True)[0].format(**row)
        outputs.append(sample)
    return outputs


def get_all_aggreact_outputs():
    outputs = []
    for idx, row in samples.reset_index().iterrows():
        sample = expand(construct_path(step="aggreact", show_replicate=False), **row)
        sampleipan = expand(
            construct_path("aggreact-ipanemap", show_replicate=False, ext=".shape"),
            **row,
        )
        outputs.append(sample)
        outputs.append(sampleipan)
    #        else:
    #            sample = expand(construct_path("qushape", ext=".qushape"), **row)
    #            outputs.append(sample)
    return outputs


# Not parallel, to be improved
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


def get_all_footprint_outputs(wildcards):
    outputs = []
    for compare in config["footprint"]["compares"]:
        outputs.extend(
            expand(
                f"{RESULTS_DIR}/{{folder}}/{{rna_id}}_footprint_{{foot_id}}.tsv",
                folder=config["folders"]["footprint"],
                rna_id=compare["rna_id"],
                foot_id=compare["id"],
                **wildcards,
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
                        f"{RESULTS_DIR}/{{folder}}/"
                        f"{{rna_id}}_pool_{{pool_id}}_{{idx}}_cond_{CONDITION}.varna",
                        folder=config["folders"]["varna"],
                        **wildcards,
                        **cond,
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
                    f"{RESULTS_DIR}/{{folder}}/"
                    f"{{rna_id}}_pool_{{pool_id}}_{{idx}}_cond_{CONDITION}.varna",
                    folder=config["folders"]["varna"],
                    rna_id=pool["rna_id"],
                    pool_id=pool["id"],
                    idx=glob.idx,
                    **cond,
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
#        sample = exppath.format(folder=folder,
#        sample=os.path.splitext(row.sample_filename)[0])
#        control = exppath.format(folder=folder,
#        sample=os.path.splitext(row.control_filename)[0])
#    return outputs
