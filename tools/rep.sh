#!/bin/sh
mkdir tmp
for i in $*
do
sed 's/2014-2025 Discrete Optimization, TU Darmstadt/2014-2026 Discrete Optimization, TU Darmstadt/g' $i > tmp/$i
done
mv tmp/* .
rmdir tmp
