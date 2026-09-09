#!/bin/bash

OUTDIR="output"
REFDIR="reference"

for name in ${OUTDIR}/*.out
do

name=${name%.out}
name=${name#${OUTDIR}/}
echo $name

cp -v $OUTDIR/$name.out  $REFDIR/$name.ref

done
