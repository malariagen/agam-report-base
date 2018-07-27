#!/bin/bash


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
export PATH=./texlive/bin/x86_64-macOSX:$PATH
export PATH=./${CONDADIR}/bin:$PATH


# use a snapshot mirror to get reproducible install
TEXREPO=https://ctanmirror.speedata.de/2017-09-01/systems/texlive/tlnet


# install miniconda
if [ ! -f miniconda.installed ]; then
    echo "[install] installing miniconda"

    # clean up any previous
    rm -rf $CONDADIR

    # download miniconda
    curl -O https://repo.continuum.io/miniconda/Miniconda3-4.3.31-MacOSX-x86_64.sh

    # install miniconda
    bash Miniconda3-4.3.31-MacOSX-x86_64.sh -b -p $CONDADIR

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
:
