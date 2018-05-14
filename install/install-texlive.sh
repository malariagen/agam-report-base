#!/bin/bash


# N.B., assume this will be executed from the root directory of a report repo.


# ensure script errors if any command fails
set -xe


# descend into dependencies directory
DEPSDIR=deps
mkdir -pv $DEPSDIR
cd $DEPSDIR


# put dependencies on the path
export PATH=./texlive/bin/x86_64-linux:$PATH


# use a snapshot mirror to get reproducible install
#TEXREPO=https://ctanmirror.speedata.de/2017-09-01/systems/texlive/tlnet
#TEXREPO=ftp://ftp.tug.org/historic/systems/texlive/2017/tlnet-final
TEXREPO=http://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2017/tlnet-final

# install texlive
if [ ! -f texlive.installed ]; then
    echo "[install] installing texlive"

    # clean up any previous
    rm -rvf texlive*
    rm -rvf install*

    # download texlive
    wget ${TEXREPO}/install-tl-unx.tar.gz

    # unpack archive
    tar zxvf install-tl-unx.tar.gz

    # run installation
    ./install-tl-*/install-tl \
        -repository=$TEXREPO \
        -profile=../agam-report-base/install/texlive.profile \
        -no-persistent-downloads \
        -no-verify-downloads

    # mark successful installation
    touch texlive.installed

else
    echo "[install] skipping texlive installation"
fi


echo "[install] installing additional texlive packages"
tlmgr option repository $TEXREPO
tlmgr_install="tlmgr install --no-persistent-downloads --no-verify-downloads --no-require-verification"
for package in $(cat ../agam-report-base/install/texlive.packages); do
    $tlmgr_install $package
done
