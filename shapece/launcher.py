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
    config: str = None, cores: int = 8, stoponerror: bool = False
):

    #launch_voila()

    if config is not None:
        config = [config]

    keepgoing = not stoponerror

    try:
        sm.snakemake(
            os.path.join("workflow", "Snakefile"),
            configfiles=config,
            cores=8,
            keepgoing=keepgoing,
            use_conda=True,
        )
    except Exception as e:
        print(e)


def main_wrapper():
    fire.Fire(launch_snakemake)
