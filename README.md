# IPANEMAP Suite

[![DOI](https://zenodo.org/badge/385251818.svg)](https://doi.org/10.5281/zenodo.14936736)


IPANEMAP Suite workflow intend to provide automation in the data treatment of SHAPE Capillary Electrophorese.

The workflow will enable you to generate structure data for an RNA fragment analysed by potentially several SHAPE experiments under different conditions (Temperature, Magnesium, Probes, etc).

## Documentation

All information about how to use this workflow can be found at :

[https://sargueil-citcom.github.io/ipasuite-docs](https://sargueil-citcom.github.io/ipasuite-docs)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

## Integration of software tools

- QuShape
- IPANEMAP
- RNAFold
- VARNA
- Custom scripts for file conversion, reactivity normalization and aggregation

## Configuration

Main parameters of `config.yaml` and their default values:

| Section | Parameter | Default | Description |
|---------|-----------|---------|-------------|
| `rawdata` | `type` | `fluo-ceq8000` | Raw data format: `fluo-ceq8000`, `fluo-fsa` or `fluo-ce` |
| `rawdata` | `path_prefix` | `resources/raw_data` | Raw files location |
| `rawdata` | `control` | `DMSO` | Control reagent name |
| `qushape` | `use_subsequence` | `false` | Concatenate data from several primers |
| `qushape` | `channels` | `RX: 0, RXS1: 2, BG: 0, BGS1: 2` | Sequencer channels (0-based) |
| `qushape` | `check_integrity` | `true` | Check the sequence stored in QuShape projects |
| `normalization` | `simple_outlier_percentile` | `98` | Values above this percentile are outliers |
| `normalization` | `simple_norm_term_avg_percentile` | `90` | Values above this percentile are averaged for normalization |
| `normalization` | `low_norm_reactivity_threshold` | `-0.3` | Values below are undetermined |
| `aggregate` | `norm_method` | `simple` | Normalization method: `simple` or `interquartile` |
| `aggregate` | `min_std` | `0.15` | Max standard deviation to accept a position directly |
| `aggregate` | `reactivity_medium` | `0.4` | Lower bound of the medium reactivity class |
| `aggregate` | `reactivity_high` | `0.7` | Lower bound of the high reactivity class |
| `ipanemap:config` | `nstructures` | `1000` | Number of sampled structures |
| `ipanemap:config` | `temperature` | `37` | Folding temperature (°C) |
| `ipanemap:config` | `slope` | `1.3` | Pseudo-energy slope (kcal/mol) |
| `ipanemap:config` | `intercept` | `-0.4` | Pseudo-energy intercept (kcal/mol) |
| `footprint:config` | `diff_thres` | `0.2` | Absolute difference threshold |
| `footprint:config` | `ratio_thres` | `0.2` | Relative difference threshold |
| `footprint:config` | `ttest_pvalue_thres` | `0.05` | t-test p-value threshold |

## Commands

| Command | Description |
|---------|-------------|
| `ipasuite init [project]` | Create a new project (`config.yaml`, `samples.tsv`, `resources/raw_data`, `results`) |
| `ipasuite prep` | Generate `samples.tsv` from raw file names (backs up the previous one) |
| `ipasuite config` | Open the configuration GUI in the browser |
| `ipasuite check` | Validate sequences, raw files, samples per run and replicate duplicates |
| `ipasuite qushape` | Create and open QuShape projects for untreated experiments |
| `ipasuite correl` | Run up to replicate aggregation and compute replicate correlations |
| `ipasuite run` | Run the full pipeline (`--dry_run`, `--rerun_incomplete`) |
| `ipasuite comparison` | Compare two dot-bracket structure models (interactive) |
| `ipasuite clean` | Remove files generated after QuShape (`--from_step <step>`) |
| `ipasuite log` | Show pipeline logs (`--clean` to erase) |
| `ipasuite unlock` | Remove the lock left by an interrupted run |
| `ipasuite convert_qushape` | Convert `.qushape` projects to `.qushapey` |

## Citation

Please cite: P. Hardouin, N. Pan, F-X. Lyonnet du Moutier, N. Chamond, Y. Ponty, S. Will, B. Sargueil, IPANEMAP Suite: a pipeline for probing-informed RNA structure modeling, NAR Genomics and Bioinformatics(2025), [Link](https://academic.oup.com/nargab/article/7/1/lqaf028/8093145?utm_source=advanceaccess&utm_campaign=nargab&utm_medium=email).
