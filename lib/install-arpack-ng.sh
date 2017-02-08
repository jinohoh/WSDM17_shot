#!/bin/sh
set -ex
git clone https://github.com/opencollab/arpack-ng.git
mkdir arpack-ng-build
cd arpack-ng-build
cmake -D BLAS_goto2_LIBRARY=../libopenblas.a ../arpack-ng
make
cd ../
ln -s arpack-ng-build/libarpack.a ./
