#!/bin/sh
mkdir tmp
for i in $*
do
sed 's/2014-2019 Discrete Optimization, TU Darmstadt/2014-2020 Discrete Optimization, TU Darmstadt/g' $i > tmp/$i
done
mv tmp/* .
rmdir tmp
