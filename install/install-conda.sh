#!/usr/bin/env bash


# N.B., assume this will be executed from the root directory of a report repo.


# ensure script errors if any command fails
set -xe


# conda setup
CONDADIR=conda
CONDANAME=${PWD##*/}

# descend into dependencies directory
DEPSDIR=deps
mkdir -pv $DEPSDIR
cd $DEPSDIR


# put dependencies on the path
export PATH=./${CONDADIR}/bin:$PATH


# install miniconda
if [ ! -f miniconda.installed ]; then
    echo "[install] installing miniconda"

    # clean up any previous
    rm -rf $CONDADIR

    if [ "$(uname)" == "Darwin" ]; then
        # Install for Mac OS X platform        
        # download miniconda
        wget --no-clobber https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

        # install miniconda
        bash Miniconda3-latest-MaxOSX-x86_64.sh -b -p $CONDADIR

    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
        # Install for GNU/Linux platform
        # download miniconda
        wget --no-clobber https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

        # install miniconda
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDADIR

    fi

    # set conda channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda update --yes conda

    # create default scientific Python environment
    conda create --yes --name=$CONDANAME python=3.6

    # mark success
    touch miniconda.installed

else
    echo "[install] skipping miniconda installation"
fi


echo "[install] installing Python packages"
source activate $CONDANAME
# ensure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# conda install
conda install --yes --file ../agam-report-base/install/conda.txt
# pypi install
pip install --no-cache-dir -r ../agam-report-base/install/pypi.txt
# clean conda caches
conda clean --yes --all
