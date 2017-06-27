#!/bin/bash

# N.B., assume this will be executed from the root directory of a report repo.

# ensure script errors if any command fails
set -eo pipefail

# descend into dependencies directory
mkdir -pv dependencies
cd dependencies

# put dependencies on the path
export PATH=./texlive/2016/bin/x86_64-linux:$PATH
export PATH=./miniconda/bin:$PATH

# install texlive
if [ ! -f texlive.installed ]; then
    echo "[install] installing texlive"

    # clean up any previous
    rm -rvf texlive*
    rm -rvf install*

    # download texlive
    wget ftp://tug.org/historic/systems/texlive/2016/tlnet-final/install-tl-unx.tar.gz

    # unpack archive
    tar zxvf install-tl-unx.tar.gz

    # run installation
    ./install-tl-20170413/install-tl \
        -repository=ftp://tug.org/historic/systems/texlive/2016/tlnet-final \
        -profile=../agam-report-base/install/texlive.profile \
        -no-persistent-downloads \
        -no-verify-downloads

    # install additional packages
    tlmgr option repository ftp://tug.org/historic/systems/texlive/2016/tlnet-final
    tlmgr install csquotes
    tlmgr install biblatex
    tlmgr install logreq
    tlmgr install xstring
    tlmgr install adjustbox
    tlmgr install collectbox
    tlmgr install todonotes
    tlmgr install siunitx
    tlmgr install tablefootnote
    tlmgr install xifthen
    tlmgr install ifmtarg
    tlmgr install preprint
    tlmgr install biber

    # mark successful installation
    touch texlive.installed

else
    echo "[install] skipping texlive installation"
fi

# install miniconda
if [ ! -f miniconda.installed ]; then
    echo "[install] installing miniconda"

    # clean up any previous
    rm -rf miniconda

    # download miniconda
    wget --no-clobber https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # install miniconda
    bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda
    conda upgrade --yes conda

    # create default scientific Python environment
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda create --yes --name=agam-report-base python=3.5

    # install Python packages
    source activate agam-report-base
    conda install --yes --file ../agam-report-base/install/conda.txt
    pip install --no-cache-dir -r ../agam-report-base/install/pypi.txt

    # clean conda caches
    conda clean --yes --all

    # mark success
    touch miniconda.installed

else
    echo "[install] skipping miniconda installation"
fi

# check to see how much space needed for cache
du -hs ./*
du -hs ./*/*
