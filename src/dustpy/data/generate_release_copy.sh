#!/bin/bash

CWDPATH="$(pwd)"
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/dusty

if [ "$1" != "" ]; then
    RELEASEDIRPATH="$1"
else
    RELEASEDIRPATH="$SCRIPTPATH/release"
fi

rm -rI $RELEASEDIRPATH/*
mkdir -p $RELEASEDIRPATH/dusty
echo '!DUSTY RELEASE' > $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/common.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/dusty.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/inout.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/kernel.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/math.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/misc.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/msg.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/nonopenmp.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/optprop.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/rdinp.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/solve_matrix.f90 >> $RELEASEDIRPATH/dusty/dusty.f90
cat $SCRIPTPATH/source/winds.f90 >> $RELEASEDIRPATH/dusty/dusty.f90

mkdir $RELEASEDIRPATH/dusty/docs
cd $SCRIPTPATH/docs
pdflatex manualV4.tex
cd $CWDPATH
cp $SCRIPTPATH/docs/manualV4.pdf $RELEASEDIRPATH/dusty/docs/manual.pdf
cp -r $SCRIPTPATH/data $RELEASEDIRPATH/dusty/data
cp $SCRIPTPATH/dusty.mas $RELEASEDIRPATH/dusty/
cp $SCRIPTPATH/userpar.inc $RELEASEDIRPATH/dusty/

echo 'all:' > $RELEASEDIRPATH/dusty/Makefile
echo -e "\t gfortran -O3 -lgomp -fopenmp -o dusty dusty.f90" >> $RELEASEDIRPATH/dusty/Makefile
echo 'gfortran -O3 -lgomp -fopenmp -o dusty.exe dusty.f90' > $RELEASEDIRPATH/dusty/compile.bat


cd $RELEASEDIRPATH
cd dusty
tar -xzf $SCRIPTPATH/examples_release.tar
cd ..
tar -cf dusty.tar dusty
cd $CWDPATH
