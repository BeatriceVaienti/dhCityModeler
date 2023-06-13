# dhCityModeller
A library to model historical cities procedurally based on GIS datasets

# Installation
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
## Creating a Virtual Environment





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