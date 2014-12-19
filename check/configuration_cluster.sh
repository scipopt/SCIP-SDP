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

# configures environment variables for cluster runs.
# to be invoked inside a check_cluster.sh script
# This script cancels the process if required variables are not correctly set

# input variables - should be passed to this script
QUEUE=$1    # the name of the cluster (M610, M620, dbg, telecom-dbg, mip-dbg, opt-low, opt)
PPN=$2      # number of cluster nodes to use
EXCLUSIVE=$3 # should cluster nodes be blocked for other users while the jobs are running?
QUEUETYPE=$4 # either 'srun' or 'qsub'

# new environment variables defined by this script:
NICE=""
ACCOUNT="mip"
CLUSTERQUEUE=$QUEUE

# check if queue has been defined
if test "$QUEUE" = ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# check if number of nodes has been defined
if test "$PPN" = ""
then
    echo Skipping test since the number of nodes has not been defined.
    exit
fi

# check whether there is enough memory on the host system, otherwise we need to submit from the target system
if test "$QUEUETYPE" = "srun"
then
    HOSTMEM=`ulimit -m`
    if test "$HOSTMEM" != "unlimited"
    then
        if [ `expr $HARDMEMLIMIT \* 1024` -gt $HOSTMEM ]
        then
            echo "Not enough memory on host system - please submit from target system (e.g. ssh opt201)."
            exit
        fi
    fi
fi

#define clusterqueue, which might not be the QUEUE, because this might be an alias for a bunch of QUEUEs
if test $CLUSTERQUEUE = "dbg"
then
    CLUSTERQUEUE="mip-dbg,telecom-dbg"
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "telecom-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "mip-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "opt-low"
then
    CLUSTERQUEUE="opt"
    NICE="--nice=10000"
elif test $CLUSTERQUEUE = "M620x"
then
    CLUSTERQUEUE="M620,M620v2,M620v3"
fi

# check if the slurm blades should be used exclusively
if test "$EXCLUSIVE" = "true"
then
    EXCLUSIVE=" --exclusive"
    if test $CLUSTERQUEUE = "opt"
    then
        CLUSTERQUEUE="M610"
    fi
else
    EXCLUSIVE=""
fi
