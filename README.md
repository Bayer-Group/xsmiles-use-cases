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

### Step 1 - Install general dependencies
Create a new virtual environment and install the dependencies defined in `requirements.txt`:

```bash
# the code has been tested with Python 3.7 
python3.7 -m venv .venv_xsmiles_usecases
source ./.venv_xsmiles_usecases/bin/activate # path to the created environment
pip3 install -r requirements.txt
```

### Step 2 - Install XSMILES and CDDD

Required: internet connection. The plugin will download RDkit MinimalLib when the JupyterLab notebook is loaded.

To be able to visualize the XAI scores, you need to install XSMILES and CDDD.

XSMILES for JupyterLab is available here: https://github.com/Bayer-Group/xsmiles-jupyterlab. You want to download the `xsmiles-0.2.2.dev0-py2.py3-none-any.whl` file from the releases page.

An unofficial package for CDDD is available in this repository: `cddd-1.2.2-py3.none.any.whl`. We packed [CDDD scripts and the CDDD default_model](https://github.com/jrwnter/cddd) into a single package to use in the notebook more easily, as well as to use with our Substitution method (`attributor.py`). Please check the `smiles_attributions` notebook to see how to we use the package and import the CDDD default model. We created this package because in certain environments, Google Drive may be blocked by firewalls.

<!-- If XSMILES cannot be installed through `pip` `requirements.txt`, check how to install it here: https://xsmiles-ipywidget. -->

<!-- ```bash
pip install xsmiles
```

or -->

```bash
pip install xsmiles-0.2.1.dev0-py2.py3-none-any.whl
pip install cddd-1.2.2-py3.none.any.whl
```

Make sure `tensorboard==1.13.1` and `tensorflow==1.13.2` were installed through `requirements.txt`, CDDD depends on them.

## Step 3 - Run JupyterLab

Run JupyterLab and choose a notebook to explore:

```bash
jupyter lab notebooks
```
