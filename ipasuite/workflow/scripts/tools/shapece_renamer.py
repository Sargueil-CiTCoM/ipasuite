#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import fire
import parse
import yaml
from ...rules import load_samples
from shutil import copy2

# from snakemake.utils import expand

# from ...rules import common
import os

def rename(
    new_file_path,
    config_file,
    full_condition_path,
    full_previous_condition_path,
    tmp_dir_path="tmp/",
    rna_id=None,
    replicate=None,
    conditions=None,
):
    with open(config_file, "r") as fd:
        config = yaml.load(fd, yaml.CLoader)
    samples = load_samples.get_samples(config)

    # wild = parse.parse(full_condition_path, new_file_path)
    print(full_condition_path)
    print(new_file_path)
    # print(wild)
    # assert wild is not None

    if "DMSO" in new_file_path:
        format_type = "control_condition"
        full_condition_path = full_condition_path.replace(
            config["format"]["condition"], config["format"]["control_condition"]
        )
        print(f"replace {full_condition_path}")
    else:
        format_type = "condition"

    condwild = parse.parse(full_condition_path, new_file_path)
    wildcards = condwild.named

    sample = load_samples.get_sample(config, samples, wildcards)
    old_file_path = full_previous_condition_path.format(
        config["format"]["previous"][format_type].format(**sample)
    )

    if not os.path.exists(tmp_dir_path):
        os.makedirs(tmp_dir_path)

    copy2(old_file_path, os.path.join(tmp_dir_path, new_file_path))


def main():
    fire.Fire(rename)


if __name__ == "__main__":
    main()
