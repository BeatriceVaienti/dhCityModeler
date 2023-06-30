# dhCityModeller
A library to model historical cities procedurally based on GIS datasets

# Installation with Conda
## Creating a Virtual Environment
the file conda_env.yml contains the environment used for development.
To install the environment run the following command:
```
conda env create --file=conda_env.yml
```
To activate the environment run the following command:
```
conda activate dhCityModeller
```
To deactivate the environment run the following command:
```
conda deactivate
```

## Enabling the versioning capabilities
To use the versioning functions it's necessary to clone in your repository the cityjson-versioning-prototipe repository as a submodule. 

to do it go in your folder, activate your virtual environment and type:
```bash
$ git submodule add https://github.com/tudelft3d/cityjson-versioning-prototype.git
```

then you need to install it. To do so move into the folder cityjson-versioning-prototype and type:

```bash
pip install --editable .
```
## Mapping the values and filling gaps
The first step of this pipeline consists in preparing your dataset to be used by the library. In particular, we need to create a geodataframe whose columns and contents respod to the provided extension. In the folder ./extension you can find the .ext.json file (containing the full extension) and a .csv file containing the flattened version of the extension, which can be employed to map the values of the dataset.
To prepare the dataset we also provide a script to fill the gaps in partially completed fields. In the notebook filling_and_mapping.ipynb and in the corresponding Colab (colab_filling_and_mapping.ipynb) you can find an example of how to use the scripts. We propose to first proceed with the filling and then with the mapping, since we acknowledge that through the mapping users may have to flatten some of the information they have to comply with the prescribed fields.

## From 2D to 3D with dhCityModeller
## If you use this library please cite the following paper:
``` 
Vaienti, B.; Petitpierre, R.; di Lenardo, I.; Kaplan, F. Machine-Learning-Enhanced Procedural Modeling for 4D Historical Cities Reconstruction. Remote Sens. 2023, 1, 0

