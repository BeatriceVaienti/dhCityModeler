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
        'geopandas @ git+https://github.com/geopandas/geopandas.git@main',
        'cadquery @ git+ https://github.com/CadQuery/cadquery.git@master',
        'ipykernel',
        'pydelatin'
    ],
)
