#!/bin/sh
mkdir tmp
for i in $*
do
sed 's/Copyright (C) 2002-2020 Zuse Institute Berlin/Copyright (C) 2002-2021 Zuse Institute Berlin/g' $i > tmp/$i
done
mv tmp/* .
rmdir tmp
