#! /usr/bin/env python3

import snakemake as sm
import os
import fire


def launch_snakemake(
    config: str = None, cores: int = 8, stoponerror: bool = False
):
    if config is not None:
        config = [config]

    keepgoing = not stoponerror

    try:
        sm.snakemake(
            os.path.join(os.path.dirname(__file__), "workflow", "Snakefile"),
            configfiles=config,
            cores=8,
            keepgoing=keepgoing,
            use_conda=True,
        )
    except Exception as e:
        print(e)


def main_wrapper():
    fire.Fire(launch_snakemake)
