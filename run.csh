#/bin/csh

cd $VFEM/conversion
mkdir -p $REFERENCE
cd $REFERENCE

mkdir -p plotXS
mkdir -p plot22
mkdir -p test

rm -rf */*

source $EMSRC/Conversion/runA.csh >& test/res.out

gzip test/*.out 
#
