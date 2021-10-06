# SHAPE-CE Snakemake workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)

Shape-CE workflow intend to provide automation in the data treatment of SHAPE Capillary Electrophorese.

It relies on :

- QuShape
- IPANEMAP
- RNAFold
- VARNA.

Custom scripts for file conversion, reactivity normalization and aggregation.

The workflow will enable you to generate structure data for a RNA fragment analysed using SHAPE with a arbitrary set of conditions (Temperature, Magnesium, Probes, etc)

## Documentation

All information about how to use this workflow can be found at :

[https://citcom-lab.github.io/shape-ce-docs](https://citcom-lab.github.io/shape-ce-docs)

## Authors

* François-Xavier Lyonnet du Moutier (@ixeft)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).



## Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.
