# Requirements 

## Operating system

The pipeline should work on any system suporting Conda. 
However, some packages (ViennaRNA, scikit-bio) are not available in Windows using conda, you will have install them manually.

Tested Operating system :
- Archlinux
- Debian 11 (Bulleyes)
- Ubuntu 16.4 to 20.4

## External software

### QuShape

In order to do the manual part of the workflow, you must install QuShape. 

The original version of QuShape relies on PythonXY, an unmaintained python2 distribution.

We provide a version that is compatible with conda, and correct few bugs with the latest python2 version.

You can find it here [https://github.com/CiTCoM-Lab/QuShape](https://github.com/CiTCoM-Lab/QuShape)

A installation tutorial is available

In addition, some part of QuShape have been converted to Python 3 and embedded with this
shape-ce pipeline.

## Dependencies 

```{note}
You don't need install manually the dependencies software. They are included 
with Miniconda and your conda environnement.
```

Main dependencies software required with pipeline :

- Miniconda
- Python >= 3.9
- Java >= 8
- Snakemake
- IPANEMAP
- ViennaRNA
- VARNA (shipped with IPANEMAP)
- Python packages: scikit-bio, pandas, fire, numpy, bsddb3
