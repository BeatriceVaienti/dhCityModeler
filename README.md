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
To prepare the dataset we also provide a script to fill the gaps in partially completed fields. In the notebook filling_and_mapping.ipynb and in the corresponding Colab (colab_filling_and_mapping.ipynb) you can find an example of how to use the scripts. We propose to first proceed with the filling and then with the mapping, since we acknowledge that through the mapping users may have to flatten some of the information they have to comply with the prescribed fields.

[Open In Colab](https://colab.research.google.com/github/BeatriceVaienti/dhCityModeler/blob/master/tests/colab_filling_and_mapping.ipynb)



## From 2D to 3D with dhCityModeller
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
ABSTRACT = {The generation of 3D models depicting cities in the past holds great potential for documentation and educational purposes. However, it is often hindered by incomplete historical data and the specialized expertise required. To address these challenges, we propose a framework for historical city reconstruction. By integrating procedural modeling techniques and machine learning models within a Geographic Information System (GIS) framework, our pipeline allows for effective management of spatial data and the generation of detailed 3D models. We developed an open-source Python module that fills gaps in 2D GIS datasets and directly generates 3D models up to LOD 2.1 from GIS files. The use of the CityJSON format ensures interoperability and accommodates the specific needs of historical models. A practical case study using footprints of the Old City of Jerusalem between 1840 and 1940 demonstrates the creation, completion, and 3D representation of the dataset, highlighting the versatility and effectiveness of our approach. This research contributes to the accessibility and accuracy of historical city models, providing tools for the generation of informative 3D models. By incorporating machine learning models and maintaining the dynamic nature of the models, we ensure the possibility of supporting ongoing updates and refinement based on newly acquired data. Our procedural modeling methodology offers a streamlined and open-source solution for historical city reconstruction, eliminating the need for additional software and increasing the usability and practicality of the process.},
DOI = {10.3390/rs15133352}
}
``` 



