function [rquant_paths envstr] = set_rquant_paths(CFG)
% rquant_paths = set_rquant_paths(CFG)
%
% -- input --
% CFG: configuration struct
%
% -- output --
% rquant_paths: paths required by rQuant
% envstr = environment variables


%%%%% genome utils %%%%%
rquant_paths = '~/svn/tools/genomes:~/svn/tools/utils:~/svn/tools/ngs:~/svn/projects/genefinding/utils:';

%%%%% rproc %%%%%
if CFG.use_rproc
  p = '~/svn/tools/rproc';
  rquant_paths = sprintf('%s:%s', p, rquant_paths);
end

%%%%% optimizer %%%%%
envstr = '';
%switch CFG.optimizer
 %case 'cplex'
  p = '/fml/ag-raetsch/share/software/matlab_tools/cplex9';
  rquant_paths = sprintf('%s:%s', p, rquant_paths);
 %case 'mosek'
  p = '/fml/ag-raetsch/share/software/mosek/6/toolbox/r2009b';
  rquant_paths = sprintf('%s:%s', p, rquant_paths);
  if CFG.use_rproc
    envstr = 'export MOSEKLM_LICENSE_FILE=/fml/ag-raetsch/share/software/mosek/6/licenses/mosek.lic; export LD_LIBRARY_PATH=/fml/ag-raetsch/share/software/mosek/6/tools/platform/linux64x86/bin';
  end
 %otherwise
  %error('unknown optimizer %s', CFG.optimizer);
%end

addpath(rquant_paths);