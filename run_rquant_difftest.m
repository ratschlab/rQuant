addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');

CFG.organism = 'elegans';
%CFG.exp = 'reads_weak_bias??' ; CFG.genes_mat_fname = 'genes_expr.mat' ;
CFG.exp = 'reads_weak_bias??_2' ; CFG.genes_mat_fname = 'genes_expr_2.mat' ;
%CFG.exp = 'genes_expr_weak_bias01_simple' ; CFG.genes_mat_fname = 'genes_expr_2.mat' ;
CFG.read_len = 75;
CFG.gene_source = 'annotation';
CFG.read_maps_select = '';

%CFG.base_dir = '/fml/ag-raetsch/nobackup/projects/rquant/simulation_2010-04-09/'
CFG.base_dir = '/fml/ag-raetsch/nobackup/projects/rquant/simulation_2010-04-05/'
CFG.read_maps_dir = sprintf('%s', CFG.base_dir);
repeat_base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/'
CFG.repeat_maps_dir = sprintf('%sannotations/%s/repeat_masker/tracks/', repeat_base_dir, CFG.organism);
CFG.repeat_maps_fn = CFG.repeat_maps_dir;

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