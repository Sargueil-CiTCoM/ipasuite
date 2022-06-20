#! /usr/bin/env python3

import snakemake as sm
import os
import fire
import subprocess

base_path = os.path.dirname(__file__)


def launch_voila(
    config: str = "config/config.yaml", samples: str = "config/samples.tsv"
):
    path = os.path.join(base_path, "configurator.ipynb")
    subprocess.run(["voila", path])


def launch_snakemake(
    action="all",
    config: str = None,
    cores: int = 8,
    stoponerror: bool = False,
    refactor: str = None,
):

    # launch_voila()

    if config is not None:
        config = [config]

    keepgoing = not stoponerror
    targets = ["all"]
    extra_config = dict()

    if refactor == "addpositions":
        extra_config["refactor_enabled"] = True
        targets = ["all_add_positions"]

    try:
        sm.snakemake(
            os.path.join(base_path, "workflow", "Snakefile"),
            configfiles=config,
            config=extra_config,
            targets=targets,
            cores=cores,
            keepgoing=keepgoing,
            use_conda=True,
            conda_prefix="~/.shapece/conda",
        )
    except Exception as e:
        print(e)


def main_wrapper():
    fire.Fire(launch_snakemake)
