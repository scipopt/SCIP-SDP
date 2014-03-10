#!/bin/sh
mkdir tmp
for i in $*
do
sed 's/Copyright (C) 2011-2012 Discrete Optimization, TU Darmstadt/Copyright (C) 2011-2014 Discrete Optimization, TU Darmstadt/g' $i > tmp/$i
done
mv tmp/* .
rmdir tmp
