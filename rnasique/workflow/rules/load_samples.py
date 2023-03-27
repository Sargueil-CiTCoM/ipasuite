from snakemake.utils import validate
import pandas as pd
import os
import warnings
import logging

def get_indexes(config, replicates_in_index=True):
    indexes = ["rna_id"]
    indexes.extend(config["conditions"])
    if replicates_in_index:
        indexes.append("replicate")
    return indexes


def get_condition_types(config):
    condition_types = {name: "str" for name in config["conditions"]}
    condition_types["id"] = "str"
    condition_types["rna_id"] = "str"
    condition_types["replicate"] = "str"
    return condition_types


def get_unindexed_samples(config):
    condition_types = get_condition_types(config)
    unindexed_samples = pd.read_csv(
        config["samples"], sep="\t", dtype=condition_types
    )  # .set_index("sample", drop=False)
    unindexed_samples = unindexed_samples.loc[
        (unindexed_samples["discard"] != "yes")
        & (unindexed_samples["discard"] is not True)
    ]
    validate(
        unindexed_samples,
        schema=f"{os.path.dirname(__file__)}/../schemas/samples.schema.yaml",
        set_default=True,
    )
    if config["qushape"]["use_subsequence"]:
        unindexed_samples["rt_begin_pos"] = (
            pd.to_numeric(unindexed_samples["rt_begin_pos"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        unindexed_samples["rt_end_pos"] = (
            pd.to_numeric(unindexed_samples["rt_end_pos"], errors="coerce")
            .fillna(0)
            .astype(int)
        )

    return unindexed_samples


def get_sample(config, samples, wildcards):
    sval = [wildcards["rna_id"]]
    sval.extend([wildcards[cond] for cond in config["conditions"]])

    if "replicate" in samples.index.names:
        sval.append(wildcards.replicate)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return samples.loc[tuple(sval)]


def get_additionnal_wildcards(wildcards):
    sample = dict(get_sample(wildcards).iloc[0])
    for w in wildcards:
        sample.pop(w, None)
    return sample


def get_samples(config, replicate_in_index=True):
    unindexed_samples = get_unindexed_samples(config)

    samples = unindexed_samples.set_index(
        get_indexes(config, replicates_in_index=replicate_in_index), drop=True
    )

    index_count = pd.Index(samples.index).value_counts()
    if replicate_in_index and len(index_count[index_count > 1]) > 0:
        logging.error("Duplicated entry in samples file")
        logging.error(f"Entries: {index_count[index_count > 1]}")
        raise Exception("Duplicated entry")
    return samples
