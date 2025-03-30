# snrnaseq-step1-preprocessing
Code used to perform dataset building and preprocessing on single-nucleus RNA sequencing data, starting from [CellBender](https://github.com/broadinstitute/CellBender/tree/master)-filtered sequencing data and [demuxlet](https://github.com/statgen/popscle) output. 

This is the first code step in a series that was used to study impacts of substance use disorders (SUDs) and HIV infection on cellular transcription in human ventral midbrain ([Wilson et al. 2025](https://doi.org/10.1101/2025.02.05.636667)).
It is intended to be used on snRNA-seq data that has been aligned (and filtered/demultiplexed as appropriate).

## Installation
This code has several dependencies; an environment with them can be built using conda (the code below was run with Anaconda 2020.11):
```
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
# build conda environment
conda create -n snrnaseq-preprocess-env python=3.8 anndata=0.8.0 \
    docutils ipython_genutils jsonschema=4.4.0 keyutils \
    matplotlib=3.7.2 numpy=1.24.4 pandas=1.4.3 python-fastjsonschema \
    pyyaml scanpy=1.9.3 scipy=1.10.1 seaborn=0.11.2 seaborn-base=0.11.2 \
    yaml scikit-misc
# activate conda environment
conda activate snrnaseq-preprocess-env
```  
              

The scripts in this repository can then be run in the command line or with a run script. An example run script and configuration file are provided; they work with input data that can
be found [here](https://doi.org/10.5281/zenodo.15109541). To use the run script once input data is present, do  
```
# make script executable
chmod +x run_preprocess.sh
run_preprocess.sh
```

## Citing this work
If you use this code in your work, please cite [Wilson et al. 2025](https://doi.org/10.1101/2025.02.05.636667):
```
@article {Wilson2025.02.05.636667,
	author = {Wilson, Alyssa M. and Jacobs, Michelle M. and Lambert, Tova Y. and Valada, Aditi and Meloni, Gregory and Gilmore, Evan and Murray, Jacinta and Morgello, Susan and Akbarian, Schahram},
	title = {Transcriptional impacts of substance use disorder and HIV on human ventral midbrain neurons and microglia},
	elocation-id = {2025.02.05.636667},
	year = {2025},
	doi = {10.1101/2025.02.05.636667},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {For people with HIV (PWH), substance use disorders (SUDs) are a prominent neurological risk factor, and the impacts of both on dopaminergic pathways are a potential point of deleterious convergence. Here, we profile, at single nucleus resolution, the substantia nigra (SN) transcriptomes of 90 postmortem donors in the context of chronic HIV and opioid/cocaine SUD, including 67 prospectively characterized PWH. We report altered microglial expression for hundreds of pro- and anti-inflammatory regulators attributable to HIV, and separately, to SUD. Stepwise, progressive microglial dysregulation, coupled to altered SN dopaminergic and GABAergic signaling, was associated with SUD/HIV dual diagnosis and further with lack of viral suppression in blood. In virologically suppressed donors, SUD comorbidity was associated with microglial transcriptional changes permissive for HIV infection. We report HIV-related downregulation of monoamine reuptake transporters specifically in dopaminergic neurons regardless of SUD status or viral load, and additional transcriptional signatures consistent with selective vulnerability of SN dopamine neurons.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2025/02/08/2025.02.05.636667},
	eprint = {https://www.biorxiv.org/content/early/2025/02/08/2025.02.05.636667.full.pdf},
	journal = {bioRxiv}
}

```






