#!/bin/sh

VERSION="2.0.0"
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
    $NAME/changelog.txt \
    $NAME/data_format.txt \
    $NAME/src/scipsdp/*.c $NAME/src/scipsdp/*.cpp $NAME/src/scipsdp/*.h \
    $NAME/src/sdpi/*.c $NAME/src/sdpi/*.cpp $NAME/src/sdpi/*.h \
    $NAME/check/*.sh \
    $NAME/check/*.awk \
    $NAME/check/testset/short.test \
    $NAME/check/testset/short.solu \
    $NAME/settings/lp_approx.set

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
