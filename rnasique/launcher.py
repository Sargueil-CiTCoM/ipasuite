#! /usr/bin/env python3

import snakemake as sm
from snakemake.utils import validate
import os
import shutil
import fire
import subprocess
from glob import glob
from ruamel.yaml import YAML
import logging
from .workflow.rules import load_samples
from .workflow.scripts.tools.utils import fasta
import filecmp
import pathlib
from colorlog import ColoredFormatter
import pandas as pd
import traceback

LOGFORMAT = "  %(log_color)s%(levelname)-8s%(reset)s | %(message)s"
DATE_FORMAT = "%b %d %H:%M:%S"
formatter = ColoredFormatter(LOGFORMAT)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
# logging.basicConfig(format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)
logger.addHandler(handler)
yaml = YAML()
base_path = os.path.dirname(__file__)

step_order = [
    "fluo-ceq8000",
    "fluo-ce",
    "subseq",
    "qushape",
    "reactivity",
    "normreact",
    "alignnormreact",
    "aggreact",
    "aggreact-ipanemap",
    "ipanemap-config",
    "ipanemap-out",
    "structure",
    "varna",
    "footprint",
]


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
                "project using `rnasique init` command or specify another path"
            )
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
    #                conda_prefix="~/.rnasique/conda",
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
        os.mkdir(project)
        shutil.copy(
            os.path.join(base_path, "config", "config.tpl.yaml"),
            os.path.join(project, "config.yaml"),
        )
        shutil.copy(
            os.path.join(base_path, "config", "samples.tpl.tsv"),
            os.path.join(project, "samples.tsv"),
        )

        os.makedirs(os.path.join(project, "resources/raw_data"))
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
                conda_prefix="~/.rnasique/conda",
            )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)

    def convert_qushape(self, rerun_incomplete=False):
        self._config = self._choose_config(self._config)
        extra_config = dict()
        extra_config["convert_qushape"] = True
        targets = ["all_convert"]

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
                force_incomplete=rerun_incomplete,
                # listrules=True,
                conda_prefix="~/.rnasique/conda",
            )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)

    def qushape(
        self,
        action="all_reactivity",
        dry_run=False,
    ):
        self._config = self._choose_config(self._config)
        targets = ["all_reactivity"]
        extra_config = dict()

        extra_config["qushape"] = {"run_qushape": True}
        self._cores = 1

        try:
            # if True:
            if self.check():
                sm.snakemake(
                    os.path.join(base_path, "workflow", "Snakefile"),
                    configfiles=[self._config] if self._config else None,
                    config=extra_config,
                    targets=targets,
                    cores=1,
                    keepgoing=self._keepgoing,
                    dryrun=dry_run,
                    use_conda=True,
                    verbose=self._verbose,
                    conda_prefix="~/.rnasique/conda",
                    rerun_triggers=["mtime"],
                )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)
            logger.error("to get more information, type : rnasique log")

    def run(
        self, action="all", dry_run=False, run_qushape=False, rerun_incomplete=False
    ):
        self._config = self._choose_config(self._config)
        targets = ["all"]
        extra_config = dict()

        if run_qushape:
            extra_config["qushape"] = {"run_qushape": True}
            self._cores = 1
        try:
            if self.check():
                # if True:
                sm.snakemake(
                    os.path.join(base_path, "workflow", "Snakefile"),
                    configfiles=[self._config] if self._config else None,
                    config=extra_config,
                    targets=targets,
                    cores=self._cores,
                    keepgoing=self._keepgoing,
                    dryrun=dry_run,
                    use_conda=True,
                    verbose=self._verbose,
                    conda_prefix="~/.rnasique/conda",
                    rerun_triggers=["mtime"],
                    force_incomplete=rerun_incomplete,
                )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)
            logger.error("to get more information, type : rnasique log")

    def log(self, step=None, clean=False, print_filename=False):
        self._config = self._choose_config(self._config)
        with open(self._config, "r") as cfd:
            config = yaml.load(cfd)
            validate(config, base_path + "/workflow/schemas/config.schema.yaml")
        if step is not None:
            pattern = f"{config['results_dir']}/logs/{step}*.log"
        else:
            pattern = f"{config['results_dir']}/logs/*.log"
        gl = glob(pattern)
        if clean:
            for file in gl:
                os.remove(file)
            if len(gl) > 0:
                logger.info("Log cleaned.")

        else:
            for file in gl:
                with open(file, "r") as fd:
                    read = fd.read()
                    if len(read) > 2:
                        if print_filename:
                            print(file)
                        print(read[:-1])

    def clean(self, from_step="reactivity", keep_log=False):
        self._config = self._choose_config(self._config)
        try:
            begin = step_order.index(from_step)
        except ValueError:
            logger.error(f"Authorized value for from_step :{step_order}")
        with open(self._config, "r") as cfd:
            config = yaml.load(cfd)
            validate(config, base_path + "/workflow/schemas/config.schema.yaml")

        for folder in step_order[begin:]:
            try:
                shutil.rmtree(
                    os.path.join(config["results_dir"], config["folders"][folder])
                )
                shutil.rmtree(
                    os.path.join(
                        config["results_dir"], "figures", config["folders"][folder]
                    )
                )
                logger.info(f"{folder} cleaned")
            except FileNotFoundError:
                logger.info(f"no folder for {folder}")

            if not keep_log:
                gl = glob(f"{config['results_dir']}/logs/{folder}*.log")
                for file in gl:
                    os.remove(file)
                if len(gl) > 0:
                    logger.info(f"{folder} logs cleaned")

    def _check_dup(self, path):
        path = pathlib.Path(path)
        if not os.path.exists(path):
            return True
        files = sorted(os.listdir(path))

        distinctFileClass = []

        for file_x in files:
            are_identical = False

            for class_ in distinctFileClass:
                are_identical = filecmp.cmp(
                    path / file_x, path / class_[0], shallow=False
                )
                if are_identical:
                    class_.append(file_x)
                    break
            if not are_identical:
                distinctFileClass.append([file_x])
        contains_dup = False
        for class_ in distinctFileClass:
            if len(class_) > 1:
                contains_dup = True
                logger.warning(
                    "files "
                    + " and ".join(class_)
                    + " are identicals. Make sur your are did not"
                    " make a mistake while preparing you replicates",
                )
        return contains_dup

    def _check_conditions(self, samples, conditions, rna_id, name, type="ipanemap"):
        sample_missing = False
        query = [f' {cname} == "{cond}" &' for cname, cond in conditions.items()] + [
            f' rna_id == "{rna_id}"'
        ]
        query = "".join(query)
        if len(samples.query(query)) == 0:
            sample_missing = True
            logger.error(
                f"{type} pool {name} cannot be handled because "
                f"no sample is available for condition: {query}"
            )
        return sample_missing

    def _check_samples(self, config, samples):
        sample_missing = False
        ipan = config["ipanemap"]["pools"]
        for pool in ipan:
            if "external_conditions" in pool:
                for cond in pool["external_conditions"]:
                    if not os.path.exists(cond["path"]):
                        logger.error(
                            f" ipanemap {cond['path']} not found in "
                            f"{pool['id']} -"
                            f" external {cond['name']}"
                        )
            for conds in pool["conditions"]:
                sample_missing = (
                    self._check_conditions(
                        samples, conds, pool["rna_id"], pool["id"], "ipanemap"
                    )
                    or sample_missing
                )
        for comp in config["footprint"]["compares"]:
            sample_missing = (
                self._check_conditions(
                    samples, comp["condition1"], comp["rna_id"], comp["id"], "footprint"
                )
                or sample_missing
            )
            sample_missing = (
                self._check_conditions(
                    samples, comp["condition2"], comp["rna_id"], comp["id"], "footprint"
                )
                or sample_missing
            )

        return sample_missing

    def check(self, verbose=False):
        seq_missing = False
        raw_missing = False
        sample_missing = False
        seq_not_fasta = False
        has_duplicate = False
        self._config = self._choose_config(self._config)
        with open(self._config, "r") as cfd:
            config = yaml.load(cfd)
            validate(config, base_path + "/workflow/schemas/config.schema.yaml")

        samples = load_samples.get_unindexed_samples(config)

        indexes = load_samples.get_indexes(config, replicates_in_index=True)
        ind_samples = samples.set_index(indexes, drop=True)
        files = (
            list(samples["probe_file"])
            + list(samples["control_file"])
            + list(samples["qushape_file"])
            + list(samples["reference_qushape_file"])
        )
        files = [
            os.path.join(config["rawdata"]["path_prefix"], str(f))
            for f in files
            if f != "" and not (f != f)
        ]

        seqs = [seq for idx, seq in config["sequences"].items()]

        for file in seqs:
            if not os.path.exists(file):
                seq_missing = True
                logger.error(f"Sequence: {file} not found")
            else:
                try:
                    fastaseqs = fasta.fasta_iter(file)
                    for fname, fseq in fastaseqs:
                        if not fseq.isupper():
                            logger.warning("Sequence contains lowercase")
                        if "T" in fseq:
                            logger.warning(
                                "Sequence contains 'T', it should contains "
                                "'U' in order to be RNA"
                            )
                except Exception as e:
                    logger.error(traceback.format_exc())
                    logger.error(e)
                    logger.error(f"Sequence: {file} does not seems to be a fasta file")
                    seq_not_fasta = True

        for file in files:
            if not os.path.exists(file):
                raw_missing = True
                logger.warning(f"Raw data : {file} not found")

        index_count = pd.Index(ind_samples.index).value_counts()
        for i, val in index_count[index_count > 1].items():
            logger.error(
                f"Found {val} Duplicated samples for"
                f" condition {dict(zip(indexes, i))}"
            )
            logger.error(
                "Every sample of the same condition must has a different replicate id. " 
                f"You must  change replicate id in the samples.tsv file"

            )
            has_duplicate = True

        self._check_dup(
            os.path.join(config["results_dir"], config["folders"]["fluo-fsa"])
        )
        self._check_dup(
            os.path.join(config["results_dir"], config["folders"]["fluo-ceq8000"])
        )

        sample_missing = self._check_samples(config, samples)
        if raw_missing or seq_missing or sample_missing or seq_not_fasta:
            logger.error("-- Problems where found when checking pipeline --")
        else:
            print("Configuration check succeed")

        return not seq_missing and not sample_missing and not seq_not_fasta and not has_duplicate


def main_wrapper():
    fire.Fire(Launcher)
