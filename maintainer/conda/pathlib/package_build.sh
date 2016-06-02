# serve as a reminder how to generate package

conda build .
conda convert $CONDA_PATH/conda-bld/linux-64/pathlib-1.0.1-py27_1.tar.bz -p linux-32
conda convert $CONDA_PATH/conda-bld/linux-64/pathlib-1.0.1-py27_1.tar.bz -p osx-64

# upload each package individually to anaconda
