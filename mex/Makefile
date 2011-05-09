include ../bin/rquant_config.sh

MEX=${MATLAB_MEX_PATH}
MKOCTFILE=${OCTAVE_MKOCT}
MATLAB_INCL=${MATLAB_INCLUDE_DIR}
SAMDIR=${SAMTOOLS_DIR}


all:	get_reads.mex get_bam_properties.mex get_reads.mexa64 get_bam_properties.mexa64
octave:	get_reads.mex get_bam_properties.mex
matlab:	get_reads.mexa64 get_bam_properties.mexa64


get_reads.mexa64:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${MEX} -g -O get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses -I$(MATLAB_INCL)

get_bam_properties.mexa64:	get_bam_properties.cpp
	rm -f *.o
	${MEX} -g -O get_bam_properties.cpp -I$(MATLAB_INCL)

get_reads.mex:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${MKOCTFILE} -g --mex get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses

get_bam_properties.mex:	get_bam_properties.cpp
	rm -f *.o
	${MKOCTFILE} -g --mex get_bam_properties.cpp

clean:	
	rm -f *.o *.mexa64 *.mex 