from setuptools import setup, find_packages

setup(
    name="geomodel_optimizer",
    version="0.5.0",
    author="Rögnvaldur L. Magnússon",
    author_email="rm@isor.is",
    description="Location optimizer for Waiwera simulations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ISOR-Geothermal/GeoModel-Optimizer",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "layermesh",
        "pywaiwera",
        "tqdm"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Hydrology"
    ],
    python_requires=">=3.10",
)
