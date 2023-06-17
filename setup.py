from setuptools import setup

setup(
    name='dhCityModeller',
    version='1.0.0',
    description='Package for GIS based procedural modelling',
    author='Beatrice Vaienti',
    packages=['modules'],  
    install_requires=[
        'numpy',
        'pandas',
        'cadquery @ git+https://github.com/dhlab-epfl/cadquery.git@master',
        'ipykernel',
        'geopandas',
        'pydelatin'
    ],
)
