%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function rquant_paths = set_rquant_paths(CFG)
% rquant_paths = set_rquant_paths(CFG)
%
% -- input --
% CFG: configuration struct
%
% -- output --
% rquant_paths: paths required by rQuant
% envstr = environment variables


%%%%% genome utils %%%%%
rquant_paths = '~/svn/tools/genomes:~/svn/tools/utils:~/svn/tools/ngs:~/svn/projects/genefinding/utils:~/svn/projects/RNASeq_galaxy/tracks:';

%%%%% rproc %%%%%
if CFG.use_rproc
  p = '~/svn/tools/rproc';
  rquant_paths = sprintf('%s:%s', p, rquant_paths);
end

addpath(rquant_paths);