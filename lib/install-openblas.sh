#!/bin/sh
set -ex
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
make
cd ../
ln -s OpenBLAS/libopenblas.a ./
#ln -s OpenBLAS/libopenblas.so ./
cd ../
