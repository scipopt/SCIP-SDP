#!/bin/sh
mkdir tmp
for i in $*
do
sed 's/Copyright (C) 2011-2017 Discrete Optimization, TU Darmstadt/Copyright (C) 2011-2018 Discrete Optimization, TU Darmstadt/g' $i > tmp/$i
done
mv tmp/* .
rmdir tmp
