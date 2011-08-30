function rquant_rproc(PAR)
% RQUANT_RPROC   Wrapper to call rQuant from rproc.
%
%   rquant_rproc(PAR)
%
%   -- input --
%   PAR contains
%     anno_dir:        directory of genes
%     track:           name of BAM file
%     output_file:     result gff3 file 
%     output_dir:      output directory
%     load_profiles:   flag that enables loading of profiles
%     profiles_fn:     name of input file that stores profiles
%     learn_profiles:  flag that enables learning of profiles
%     profiles_fn_out: name of output file that stores profiles
%     CFG:             configuration struct
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%


%%%%% paths %%%%%
% rQuant paths
global RQUANT_PATH RQUANT_SRC_PATH
% rQuant version
global RQUANT_VERSION
% SAMTools path
global SAMTOOLS_DIR
RQUANT_VERSION='2.1';
RQUANT_PATH='/fml/ag-raetsch/home/bohnert/svn/releases/rQuant/trunk';
RQUANT_SRC_PATH='/fml/ag-raetsch/home/bohnert/svn/projects/rquant/src-dev';
SAMTOOLS_DIR='/fml/ag-raetsch/share/software/samtools';
addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');
%RQUANT_VERSION = getenv('RQUANT_VERSION');
%RQUANT_PATH = getenv('RQUANT_PATH');
%RQUANT_SRC_PATH = getenv('RQUANT_SRC_PATH');
%SAMTOOLS_DIR = getenv('SAMTOOLS_DIR');
if isempty(RQUANT_VERSION) || isempty(RQUANT_PATH) || isempty(RQUANT_SRC_PATH) || isempty(SAMTOOLS_DIR)
  error('rQuant paths are not set. At least one of these variables is empty: RQUANT_VERSION, RQUANT_PATH, RQUANT_SRC_PATH, SAMTOOLS_DIR.');
end
rquant(PAR.anno_dir, PAR.track, PAR.output_file, PAR.output_dir, PAR.load_profiles, PAR.profiles_fn, PAR.learn_profiles, PAR.profiles_fn_out, PAR.CFG);