# Installation

- Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

- Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).

- Clone your repository 
```bash
git clone --recurse-submodules https://path.to/your/repository/shape-ce

```

- From repository root, create conda environnement
```bash
conda env create --name snakemake --file workflow/envs/snakemake.yml
```
- Start conda environement
```
source activate snakemake
# or 
conda activate snakemake
```
