# Requirements 

## Operating system

The pipeline should work on any system suporting Conda. 
However, some packages (ViennaRNA, scikit-bio) are not available in Windows using conda, you will have install them manually.

Tested Operating system :
- Archlinux
- Debian 11 (Bulleyes)
- Ubuntu 16.4 to 20.4


## Software

```{note}
You don't need install manually all the requirements software. They are included 
with Miniconda and your conda environnement.
```

All software requirements are installed through conda and the bioconda and conda-forge
channels

Here is an list of main software required with pipeline :

- Miniconda
- Python >= 3.9
- Java >= 8
- Snakemake
- QuShape ([https://github.com/CiTCoM-Lab/QuShape ](https://github.com/CiTCoM-Lab/QuShape) this fork correct some small bugs and provide a conda environnement for it)
- IPANEMAP
- ViennaRNA
- VARNA (shipped with IPANEMAP)
- Python packages: scikit-bio, pandas, fire, numpy, bsddb3
