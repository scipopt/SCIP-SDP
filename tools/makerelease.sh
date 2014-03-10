#!/bin/sh

VERSION="1.0"
NAME="scipsdp-$VERSION"
rm -f $NAME

# add link in order to create right names
ln -s . $NAME

# run git status to clean the dirty git hash
git status
echo ""

tar --no-recursion --ignore-failed-read -cvzhf $NAME.tgz \
    --exclude="*~" \
    --exclude=".*" \
    $NAME/Makefile \
    $NAME/doc/scipsdp.dxy \
    $NAME/doc/xternal.c \
    $NAME/src/*.c $NAME/src/*.cpp $NAME/src/*.h

echo ""
echo "check version numbers in doc/scipsdp.dxy, doc/xternal.c, and Makefile:"
echo -n "scipsdp.dxy: "
grep "^PROJECT\_NUMBER" doc/scipsdp.dxy
echo -n "xternal.c:   "
grep "\@version" doc/xternal.c
echo -n "Makefile:    "
grep "^SCIPSDPVERSION" Makefile
echo ""

# remove link again
rm -f $NAME
