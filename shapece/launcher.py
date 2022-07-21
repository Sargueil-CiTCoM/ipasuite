#! /usr/bin/env python3

import snakemake as sm
import os
import shutil
import fire
import subprocess

base_path = os.path.dirname(__file__)


class Launcher(object):
    def __init__(
        self,
        config: str = None,
        cores: int = 8,
        stoponerror: bool = False,
        verbose: bool = False,
    ):
        self._config = config
        self._cores = cores
        self._keepgoing = not stoponerror
        self._verbose = verbose

    def _choose_config(self, config):
        config = (
            config
            if config is not None
            else (
                "config.yaml" if os.path.exists("config.yaml") else "config/config.yaml"
            )
        )
        if not os.path.exists(config):
            raise fire.core.FireError(
                f"{config} file does not exist, please init your "
                "project using `shapece init` command or specify another path"
            )
        print(config)
        return config

    def config(self, dev=False):
        self._config = self._choose_config(self._config)
        path = os.path.join(base_path, "configurator.ipynb")
        env = os.environ.copy()
        env["CONFIG_FILE_PATH"] = os.path.join(self._config)
        env["PROJECT_PATH"] = os.path.join(os.getcwd())
        if dev:
            subprocess.Popen(["jupyter-notebook", path], env=env)
        else:
            subprocess.Popen(["voila", path], env=env)

    #    def report(self):
    #        targets = ["all"]
    #        extra_config = dict()
    #        report_path =os.path.join(os.getcwd(),"report.html")
    #        try:
    #            sm.snakemake(
    #                os.path.join(base_path, "workflow", "Snakefile"),
    #                configfiles=[self._config],
    #                config=extra_config,
    #                targets=targets,
    #                cores=self._cores,
    #                report=report_path,
    #                keepgoing=self._keepgoing,
    #                use_conda=True,
    #                verbose=self._verbose,
    #                conda_prefix="~/.shapece/conda",
    #            )
    #        except Exception as e:
    #            print(e)
    #        webbrowser.open_new_tab(report_path)
    def report(self, dev=False):
        self._config = self._choose_config(self._config)
        path = os.path.join(base_path, "report.ipynb")
        env = os.environ.copy()
        env["CONFIG_FILE_PATH"] = os.path.join(self._config)
        env["PROJECT_PATH"] = os.path.join(os.getcwd())
        if dev:
            subprocess.Popen(["jupyter-notebook", path], env=env)
        else:
            subprocess.Popen(["voila", path], env=env)

    def init(self, project: str):

        if os.path.exists(project):
            fire.core.FireError(f"{project} folder already exists")

        os.mkdir(os.path.join(project, "config"))
        shutil.copy(
            os.path.join(base_path, "config", "config.tpl.yaml"),
            os.path.join(project, "config.yaml"),
        )
        shutil.copy(
            os.path.join(base_path, "config", "samples.tpl.tsv"),
            os.path.join(project, "samples.tsv"),
        )

        os.makedirs(os.path.join(project, "resources"))
        os.mkdir(os.path.join(project, "results"))

    def refactor(self, action: str = "addpositions"):
        extra_config = dict()

        if action == "addpositions":
            extra_config["refactor_addpositions"] = True
            targets = ["all_add_positions"]
        elif action == "rename":
            extra_config["refactor_rename"] = True
            targets = ["all_rename"]

        else:
            raise fire.core.FireError(f"invalid refactor option {action}")

        try:
            sm.snakemake(
                os.path.join(base_path, "workflow", "Snakefile"),
                configfiles=[self._config] if self._config else None,
                config=extra_config,
                targets=targets,
                cores=self._cores,
                keepgoing=self._keepgoing,
                use_conda=True,
                verbose=self._verbose,
                # listrules=True,
                conda_prefix="~/.shapece/conda",
            )
        except Exception as e:
            print(e)

    def run(
        self,
        action="all",
    ):
        targets = ["all"]
        extra_config = dict()

        try:
        #if True:
            sm.snakemake(
                os.path.join(base_path, "workflow", "Snakefile"),
                configfiles=[self._config] if self._config else None,
                config=extra_config,
                targets=targets,
                cores=self._cores,
                keepgoing=self._keepgoing,
                use_conda=True,
                verbose=self._verbose,
                conda_prefix="~/.shapece/conda",
            )
        except Exception as e:
            print(e)


def main_wrapper():
    fire.Fire(Launcher)
