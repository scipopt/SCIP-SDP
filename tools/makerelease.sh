#!/bin/sh

VERSION="3.1.2"
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
    $NAME/doc/layout.xml \
    $NAME/instances/* \
    $NAME/INSTALL \
    $NAME/license.txt \
    $NAME/changelog.txt \
    $NAME/sdpa_format.txt \
    $NAME/src/scipsdp/*.c $NAME/src/scipsdp/*.cpp $NAME/src/scipsdp/*.h \
    $NAME/src/sdpi/*.c $NAME/src/sdpi/*.cpp $NAME/src/sdpi/*.h \
    $NAME/src/scipsdpgithash.c \
    $NAME/check/testset/short.test \
    $NAME/check/testset/short.solu \
    $NAME/settings/lp_approx.set \
    $NAME/CMakeLists.txt \
    $NAME/src/CMakeLists.txt \
    $NAME/cmake/Modules/*

echo ""
echo "check version numbers in doc/scipsdp.dxy, doc/xternal.c, Makefile, and CMakeLists.txt:"
echo -n "scipsdp.dxy: "
grep "^PROJECT\_NUMBER" doc/scipsdp.dxy
echo -n "xternal.c:   "
grep "\@version" doc/xternal.c
echo -n "Makefile:    "
grep "^SCIPSDPVERSION" Makefile
echo ""

# remove link again
rm -f $NAME
