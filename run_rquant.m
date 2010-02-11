addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');

CFG.organism = 'elegans';
CFG.exp = 'nGASP-Train'; % 'elegans_pair_15'
CFG.read_len = 36;
CFG.gene_source = 'annotation';
CFG.read_maps_select = '_el8_mm1'; % ''

CFG

rproc_memreq                = 5000;
rproc_par.priority          = 500;
rproc_par.express           = 0;
rproc_par.immediately_bg    = 0;
rproc_par.immediately       = 0;
rproc_par.arch              = 64;
rproc_par.identifier        = '';
rproc_par.verbosity         = 0;
rproc_time                  = 72*60; % mins
if 0
  rproc_par.envstr = 'export MOSEKLM_LICENSE_FILE=/fml/ag-raetsch/home/bohnert/tmp/mosek.lic; export LD_LIBRARY_PATH=/fml/ag-raetsch/share/software/mosek/5/tools/platform/linux64x86/bin';
  rproc_par.identifier = sprintf('rq.%s-', organism(1:2));
  fprintf(1, 'Submitting job %i (%s) to cluster\n', cnt, rproc_par.identifier);
  job(cnt) = rproc('rquant', CFG, rproc_memreq, rproc_par, rproc_time);
else
  rquant(CFG);
end