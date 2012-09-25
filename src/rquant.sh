#!/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
#

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/rquant_config.sh

echo
echo ${PROG}: This program is part of the rQuant version $RQUANT_VERSION.
echo
echo rQuant determines the abundance of multiple transcripts per 
echo gene locus from RNA-Seq measurements.
echo 

if [ -z "$10" -o "$1" == '--help' ];
then
  echo Usage: $0
  echo "   or:" $0 --help
  false
fi 

ANNO_INPUT=${1}
ANNO_FORMAT=${2}
GENES_FN=${3}
BAM_INPUT=${4}
RQUANT_RES_FILE=${5}
RQUANT_RES_DIR=${6}
LOAD_PROFILES=${7}
PROFILES_FN=${8}
LEARN_PROFILES=${9}
PROFILES_FN_OUT=${10}

mkdir -p $RQUANT_RES_DIR


echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

if [ "$ANNO_FORMAT" == '0' ]
then
    echo load the genome annotation in GFF3 format and create an annotation object
    echo
    if [[ ! -f ${GENES_FN} || ! -s ${GENES_FN} ]]
    then
	export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
	${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} all all ${GENES_FN}
	if [[ ! -s ${GENES_FN} && -f ${GENES_FN}.mat ]]
	then
	    mv ${GENES_FN}.mat ${GENES_FN}
	fi
    fi    
fi

if [ "$ANNO_FORMAT" == '1' ]
then
    echo load the genome annotation in AGS format
    echo
    GENES_FN=${ANNO_INPUT}
fi
${DIR}/../bin/genes_cell2struct ${GENES_FN}
ln -s ${GENES_FN} ${RQUANT_RES_DIR}/genes.mat


echo
echo %%%%%%%%%%%%%%%%%%%%%
echo % 2. Quantification %
echo %%%%%%%%%%%%%%%%%%%%%
echo

echo quantify transcripts using given alignments
echo 
${DIR}/../bin/rquant ${RQUANT_RES_DIR} ${BAM_INPUT} ${RQUANT_RES_FILE} ${RQUANT_RES_DIR}/ ${LOAD_PROFILES} ${PROFILES_FN} ${LEARN_PROFILES} ${PROFILES_FN_OUT}


echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
