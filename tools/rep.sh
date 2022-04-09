#!/bin/sh
mkdir tmp
for i in $*
do
sed 's/2014-2021 Discrete Optimization, TU Darmstadt/2014-2022 Discrete Optimization, TU Darmstadt/g' $i > tmp/$i
done
mv tmp/* .
rmdir tmp
