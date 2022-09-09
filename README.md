# XSMILES - JupyterLab example notebooks

Examples of pipelines from *models & explanations* to *visualizations*.

Available notebooks in `notebooks/`:

 - Visualizing Gasteiger Charges (Simple example): `notebooks/atom_attributions_gasteiger_charges.ipynb`

- Comparing LogD and bioconcetration factor attributions (Loading attributions from JSON): `TBD`.

- Comparing LogP attributions from different methods (from ML Models to Attributions and Visualization): `notebooks/smiles_attributions_for_logp.ipynb`

## Please Cite

If you use XSMILES, the use cases, its code, or the generated explanations, please cite our article:

**Article in preparation**, a preprint should appear mid September!

```prose
Heberle, H., Zhao, L., Schmidt, S., Wolf, T., & Heinrich, J. (2022). XSMILES: interactive visualization for molecules, SMILES and XAI scores. Article in preparation.
```

```BibTeX
@article{Heberle2022XSMILES,
author = {Heberle, Henry and Zhao, Linlin and Schmidt, Sebastian and Wolf, Thomas and Heinrich, Julian},
doi = {},
journal = {Article in preparation},
month = {},
number = {},
pages = {},
title = {{XSMILES: interactive visualization for molecules, SMILES and XAI scores}},
volume = {},
year = {2022}
}
```

![JupyterLab Notebook](/screenshot.png)

## XSMILES for Javascript, KNIME, and How to use it

- [XSMILES main project](https://github.com/Bayer-Group/xsmiles)

- [XSMILES for JupyterLab](https://github.com/Bayer-Group/xsmiles-jupyterlab)

## How to run the notebook

### Step 1 - Install general dependencies and XSMILES
Create a new virtual environment and install the dependencies defined in `requirements.txt`:

```bash
# the code has been tested with Python 3.7, it's a dependency from CDDD
python3.7 -m venv .venv_xsmiles_usecases
source ./.venv_xsmiles_usecases/bin/activate # path to the created environment
pip3 install -r requirements.txt
```

### Step 2 - Install CDDD 

An unofficial package for CDDD is available in this repository: `cddd-1.2.2-py3.none.any.whl`. We packed [CDDD scripts and the CDDD default_model](https://github.com/jrwnter/cddd) into a single package to use in the notebook more easily, as well as to use with our Substitution method (`attributor.py`). Please check the `smiles_attributions` notebook to see how to we use the package and import the CDDD default model. We created this package because in certain environments, Google Drive may be blocked by firewalls.

```bash
pip install cddd-1.2.2-py3.none.any.whl
```

Make sure `tensorboard==1.13.1` and `tensorflow==1.13.2` were installed correctly through `requirements.txt`, CDDD depends on them, as well as on `python <= 3.7`.

You can use XSMILES for JupyterLab with newer versions of python. This dependency on Python 3.7 is here only for the CDDD model to work.

## Step 3 - Run JupyterLab

Run JupyterLab and choose a notebook to explore:

```bash
jupyter lab notebooks
```

## Notes

### XSMILES from .whl file

If you don't want to install XSMILES from pipy (requirements.txt), you can install the .whl file available [here](https://github.com/Bayer-Group/xsmiles-jupyterlab)

```bash
pip install xsmiles-0.2.1.dev0-py2.py3-none-any.whl
```

### Internet connection is a requirement

The plugin will download RDkit MinimalLib when the JupyterLab notebook is loaded.


