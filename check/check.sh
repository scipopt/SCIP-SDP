#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
TSTNAME=$1
BINNAME=$2
SETNAMES=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
DISPFREQ=${10}
CONTINUE=${11}
LOCK=${12}
VERSION=${13}
LPS=${14}
VALGRIND=${15}
CLIENTTMPDIR=${16}
OPTCOMMAND=${17}

# check if all variables defined (by checking the last one)
if test -z $OPTCOMMAND
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAMES       = $SETNAMES"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "LOCK          = $LOCK"
    echo "VERSION       = $VERSION"
    echo "LPS           = $LPS"
    echo "VALGRIND      = $VALGRIND"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "OPTCOMMAND    = $OPTCOMMAND"
    exit 1;
fi

# call routines for creating the result directory, checking for existence
# of passed settings, etc
TIMEFORMAT="sec"
MEMFORMAT="kB"
. ./configuration_set.sh $BINNAME $TSTNAME $SETNAMES $TIMELIMIT $TIMEFORMAT $MEMLIMIT $MEMFORMAT $VALGRIND

INIT="true"
COUNT=0
for INSTANCE in `cat testset/$TSTNAME.test` DONE
do
    COUNT=`expr $COUNT + 1`

    # loop over settings
    for SETNAME in ${SETTINGSLIST[@]}
    do
    # infer the names of all involved files from the arguments
        p=0 # currently, noone uses permutations here
        PERMUTE=0
        QUEUE=`hostname`
        # infer the names of all involved files from the arguments
        . ./configuration_logfiles.sh $INIT $COUNT $INSTANCE $BINID $PERMUTE $SETNAME $TSTNAME $CONTINUE $QUEUE  $p

        if test "$INSTANCE" = "DONE"
        then
            #echo $EVALFILE
            ./evalcheck_cluster.sh -r $EVALFILE
            continue
        fi
        # check if problem instance exists
        SCIP_INSTANCEPATH=$SCIPPATH
        for IPATH in ${POSSIBLEPATHS[@]}
        do
            echo $IPATH
            if test "$IPATH" = "DONE"
            then
                echo "input file $INSTANCE not found!"
                SKIPINSTANCE="true"
            elif test -f $IPATH/$INSTANCE
            then
                SCIP_INSTANCEPATH=$IPATH
                break
            fi

        done

        if test "$SKIPINSTANCE" = "true"
        then
            continue
        fi

        INSTANCE=${SCIP_INSTANCEPATH}/${INSTANCE}

        # overwrite the tmp file now
        # call tmp file configuration for SCIP
        . ./configuration_tmpfile_setup_scip.sh $INSTANCE $SCIPPATH $TMPFILE $SETNAME $SETFILE $THREADS $SETCUTOFF $FEASTOL $TIMELIMIT $MEMLIMIT $NODELIMIT $LPS $DISPFREQ $OPTCOMMAND $SOLUFILE

        # additional environment variables needed by run.sh
        export SOLVERPATH=$SCIPPATH
        export EXECNAME=${VALGRINDCMD}$SCIPPATH/../$BINNAME
        export BASENAME=$FILENAME
        export FILENAME=$INSTANCE
        export CLIENTTMPDIR
        echo Solving instance $INSTANCE with settings $SETNAME, hard time $HARDTIMELIMIT, hard mem $HARDMEMLIMIT
        bash -c "ulimit -t $HARDTIMELIMIT s; ulimit -v $HARDMEMLIMIT k; ulimit -f 200000; ./run.sh"
        #./run.sh
    done
    INIT="false"
done
