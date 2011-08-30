#/bin/bash

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

. ../bin/rquant_config.sh

PROG=`basename $0`

echo
echo ${PROG}: This program is part of the rQuant version $RQUANT_VERSION.
echo
echo rQuant determines the abundance of multiple transcripts per 
echo gene locus from RNA-Seq measurements.
echo 

if [ -z "$1" -o "$1" == '--help' ];
then
  echo Usage: $0 small\|big
  echo "   or:" $0 --help
  false
fi 
if [ "$1" != 'small' -a "$1" != 'big' ];
then
  echo invalid parameter
  false
fi

if [ "$1" == 'small' ];
then
  echo Note: Running this script takes about 1 minute \(on a single CPU\).
  FASTA_INPUT=data/nGASP-Train-I.fasta  
  GFF3_INPUT=data/nGASP-Train-I.gff3
  SAM_INPUT=data/nGASP-Train-I.sam
  BAM_INPUT=data/nGASP-Train-I.bam
  EXP=nGASP-Train-I
fi
if [ "$1" == 'big' ];
then
  echo Note: Running this script takes about 5 minutes \(on a single CPU\).
  FASTA_INPUT=data/nGASP-Train.fasta 
  GFF3_INPUT=data/nGASP-Train.gff3
  SAM_INPUT=data/nGASP-Train.sam
  BAM_INPUT=data/nGASP-Train.bam
  EXP=nGASP-Train
fi
LOAD_PROFILES="0"
PROFILES_FN="dummy_pr_in"
LEARN_PROFILES="0"
PROFILES_FN_OUT="dummy_pr_out"

RESULTDIR=./results-$1
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

GENES_RES_DIR=$RESULTDIR/elegans.anno
mkdir -p $GENES_RES_DIR
GENES_FN=${GENES_RES_DIR}/genes.mat

echo 1a. load the genome annotation in GFF3 format \(format version = wormbase\), create an annotation object #\[Log file in ${GENES_RES_DIR}/elegans-gff2anno.log\]
if [ ! -f ${GENES_FN} ]
then
    export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
    ${PYTHON_PATH} -W ignore::FutureWarning ../tools/ParseGFF.py ${GFF3_INPUT} all all ${GENES_FN} #> ${RESULTDIR}/elegans-gff2anno.log
    ../bin/genes_cell2struct ${GENES_FN}
fi

echo 1b. convert the alignments in SAM format to BAM format
../tools/./sam_to_bam.sh ${SAMTOOLS_DIR} data ${FASTA_INPUT} ${SAM_INPUT}

echo
echo %%%%%%%%%%%%%%%%%%%%%
echo % 2. Quantification %
echo %%%%%%%%%%%%%%%%%%%%%
echo

RQUANT_RES_DIR=$RESULTDIR/rquant
mkdir -p $RQUANT_RES_DIR

echo quantify transcripts using given alignments \(log file in ${RESULTDIR}/elegans-rquant.log\)
../bin/rquant ${GENES_RES_DIR} ${BAM_INPUT} ${RQUANT_RES_DIR}/${EXP}_rquant.gff3 ${RQUANT_RES_DIR}/ ${LOAD_PROFILES} ${PROFILES_FN} ${LEARN_PROFILES} ${PROFILES_FN_OUT} > ${RESULTDIR}/elegans-rquant.log

echo
echo Quantification result can be found in ${RQUANT_RES_DIR}/${EXP}_rquant.gff3
echo

echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
