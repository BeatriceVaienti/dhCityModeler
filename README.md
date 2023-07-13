# dhCityModeler
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
To prepare the dataset we also provide a script to fill the gaps in partially completed fields. In the notebook /tests/filling_and_mapping.ipynb and in the corresponding Colab (colab_filling_and_mapping.ipynb) you can find an example of how to use the scripts. We propose to first proceed with the filling and then with the mapping, since we acknowledge that through the mapping users may have to flatten some of the information they have to comply with the prescribed fields.

[Open In Colab](https://colab.research.google.com/github/BeatriceVaienti/dhCityModeler/blob/master/tests/colab_filling_and_mapping.ipynb)


## From 2D to 3D with dhCityModeller
Once you have a geojson that has been mapped to the prescribed fields you're all set to proceed with the encoding in the CityJSON format and the 3D generation of geometry (NB: it's not mandatory to fill every field! If some necessary field it's missing it will be automatically filled and marked as such in its paradata).
The core principle of our library is that 3D geometry is the result of the combination between a footprint and a set of parameters that describe the geometry. The geometry generation is built using [CadQuery](https://github.com/CadQuery/cadquery), a library for 3D CAD modeling.
You can test the encoding and modeling functions by using the provided notebook (/tests/3D_generation.ipynb) or the corresponding Google Colab:

[Open In Colab](https://colab.research.google.com/github/BeatriceVaienti/dhCityModeler/blob/master/tests/colab_3D_generation.ipynb)) 

## If you use this library please cite the following paper:
``` 
@Article{rs15133352,
AUTHOR = {Vaienti, Beatrice and Petitpierre, Rémi and di Lenardo, Isabella and Kaplan, Frédéric},
TITLE = {Machine-Learning-Enhanced Procedural Modeling for 4D Historical Cities Reconstruction},
JOURNAL = {Remote Sensing},
VOLUME = {15},
YEAR = {2023},
NUMBER = {13},
ARTICLE-NUMBER = {3352},
URL = {https://www.mdpi.com/2072-4292/15/13/3352},
ISSN = {2072-4292},
DOI = {10.3390/rs15133352}
}
``` 

## The Historical CityJSON Extension
In the folder /extension you can find the extension file for CityJSON and its flattened form that we employ as a tabular form in the GeoJSON files that we want to encode.
This constitutes a work in progress that we are constantly updating to accomodate better the needs of our historical models and their interoperability with existing schema. You can find the principles behind it in the following paper, even though the version in this repository represents the latest version of it:

``` 
@inproceedings{10.1145/3557919.3565813,
author = {Vaienti, Beatrice and Guhennec, Paul and di Lenardo, Isabella},
title = {A Data Structure for Scientific Models of Historical Cities: Extending the CityJSON Format},
year = {2022}, isbn = {9781450395335},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3557919.3565813},
doi = {10.1145/3557919.3565813},
booktitle = {Proceedings of the 6th ACM SIGSPATIAL International Workshop on Geospatial Humanities},
pages = {20–23},
numpages = {4}, keywords = {uncertainty visualisation, version controlling, digitisation, architectural reconstruction, 4D cities, historical validation, historical modelling},
location = {Seattle, Washington}, series = {GeoHumanities '22} }
``` 

