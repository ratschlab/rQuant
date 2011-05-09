function rquant_config
% RQUANT_CONFIG   Sets a few global variables with system dependent paths.
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%


% rQuant paths
global RQUANT_PATH RQUANT_SRC_PATH RQUANT_TMP_PATH RQUANT_LD_LIBRARY_PATH

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

% optimizer related paths
global OPTIMIZER OPTIMIZER_PATH OPTIMIZER_TOOLBOX_PATH

% SAMTools path
global SAMTOOLS_DIR

% configuration (adapt to the user's configuration)
RQUANT_PATH = getenv('RQUANT_PATH');
RQUANT_SRC_PATH = getenv('RQUANT_SRC_PATH');
RQUANT_TMP_PATH = getenv('RQUANT_TMP_PATH');
RQUANT_LD_LIBRARY_PATH = getenv('RQUANT_LD_LIBRARY_PATH');
INTERPRETER = getenv('INTERPRETER');
MATLAB_BIN_PATH = getenv('MATLAB_BIN_PATH');
OCTAVE_BIN_PATH = getenv('OCTAVE_BIN_PATH');
OPTIMIZER = getenv('OPTIMIZER');
OPTIMIZER_PATH = getenv('OPTIMIZER_PATH');
OPTIMIZER_TOOLBOX_PATH = getenv('OPTIMIZER_TOOLBOX_PATH');
SAMTOOLS_DIR = getenv('SAMTOOLS_DIR');

% switch off a few expected warnings
addpath(sprintf('%s/tools', RQUANT_PATH));
engine = determine_engine();
if isequal(engine, 'octave'),
  warning('off', 'mgene:no_model');
  warning('off', 'Octave:precedence-change');
  warning('off', 'Octave:function-name-clash');
  warning('off', '');
  warning('off', 'Octave:num-to-str');
  warning('off', 'Octave:function-name-clash');
  warning('off', 'Octave:divide-by-zero');
  warning('off', 'Octave:future-time-stamp');
  warning('off', 'load_genomic:contig_boundary');
  warning('off', 'merge_blocks:multiple_truths');
  warning('off', 'solve_qp:constraints');
  warning('off', 'correct_cleave_and_tss_pos:find_nearest_pos');
  warning('off', 'Octave:assign-as-truth-value');
else
  warning('off', 'MATLAB:typeaheadBufferOverflow');
end

% make sure no process stops with a debug prompt
global g_ignore_keyboard
g_ignore_keyboard = 1;
