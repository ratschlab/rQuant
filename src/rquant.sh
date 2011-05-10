#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2010 Max Planck Society
#

set -e 

PROG=`basename $0`
DIR=`dirname $0`

echo
echo ${PROG}: This program is part of the rQuant version 1.0.
echo
echo rQuant determines the abundance of multiple transcripts per 
echo gene locus from RNA-Seq measurements.
echo 

if [ -z "$11" -o "$1" == '--help' ];
then
  echo Usage: $0
  echo "   or:" $0 --help
  false
fi 

GFF3_INPUT=${1}
BAM_INPUT=${2}
RQUANT_RES_FILE=${3}
RQUANT_RES_DIR=${4}
LOAD_PROFILES=${5}
PROFILES_FN=${6}
INTRON_DISTS_FN=${7}
LEARN_PROFILES=${8}
NUM_ITER=${9}
PROFILES_FN_OUT=${10}
INTRON_DISTS_FN_OUT=${11}
GENES_FN=${RQUANT_RES_DIR}/genes.mat

mkdir -p $RQUANT_RES_DIR

. ${DIR}/../bin/rquant_config.sh

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
 
echo load the genome annotation in GFF3 format and create an annotation object
if [ ! -f ${GENES_FN} ]
then
    export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
    ${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${GFF3_INPUT} all all ${GENES_FN}
fi

echo
echo %%%%%%%%%%%%%%%%%%%%%
echo % 2. Quantification %
echo %%%%%%%%%%%%%%%%%%%%%
echo

echo quantify transcripts using given alignments
${DIR}/../bin/rquant ${RQUANT_RES_DIR} ${BAM_INPUT} ${RQUANT_RES_FILE} ${RQUANT_RES_DIR}/ ${LOAD_PROFILES} ${PROFILES_FN} ${INTRON_DISTS_FN} ${LEARN_PROFILES} ${NUM_ITER} ${PROFILES_FN_OUT} ${INTRON_DISTS_FN_OUT}

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
