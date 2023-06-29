from setuptools import setup

setup(
    name='dhCityModeler',
    version='1.0.0',
    description='Package for GIS based procedural modelling and ML based data completion.',
    author='Beatrice Vaienti',
    packages=['modules'],  
    install_requires=[
        'numpy',
        'pandas',
        'geopandas @ git+https://github.com/geopandas/geopandas.git@main',
        'cadquery-ocp',
        'cadquery @ git+https://github.com/CadQuery/cadquery.git@master',
        'ipykernel',
        'pydelatin',
        'scikit-learn>=1.2'
    ],
)
