#!/bin/bash


cd gsl-2.6
sudo chmod 777 *
./configure
make
make install
cd ..
apt-get install -y libboost-all-dev
Rscript install_packages.R
