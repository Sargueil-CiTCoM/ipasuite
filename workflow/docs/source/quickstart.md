# Quickstart 
In order to use the pipeline you will have to follow those steps :

1. Install conda and set up your git repository
2. Choose your conditions and pipeline parameters in `config/config.yaml`
3. Fill up your experiments in `config/samples.tsv` specifing where to find your raw data
   and QuShape projects if you already treated some experiments
4. Launch `./shape-ce.sh`
5. Open your generated .qushape file in results/2-qushape and perform the qushape
   analysis and save the complete analysis
6. Launch again `./shape-ce.sh`

You will find your structure data in results/5.3-structure and results/05.4-varna treated
with ipanemap as set up in your `config/config.yaml`

Data will be organized using a `config/samples.tsv` and configured using `config/config.yaml` , and be feeded to IPANEMAP with a fully configured set of conditions.
