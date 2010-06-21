addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% directories from which to load read data and genes %%%%%
CFG.organism = 'elegans';
CFG.exp = 'C_elegans_L3';
CFG.read_len = 75;
CFG.gene_source = 'annotation';

CFG.tracks_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/read_maps/elegans/C_elegans_L3/';
CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/repeat_masker/tracks/';

%%%%% genes %%%%% 
CFG.gene_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/c_elegans.WS200.new.gann/all_all_genes_merged_strands.mat';

%%%%% genome info %%%%%
CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/elegans/elegans.gio/genome.config');
% contig length
for c = 1:length(CFG.genome_info.flat_fnames),
  d = dir(CFG.genome_info.flat_fnames{c});
  CFG.chr_len(c) = d(1).bytes;
end

%%%%% alignments %%%%%
for c = 1:length(CFG.genome_info.flat_fnames),
  CFG.tracks_fn{c} = {sprintf('%s%s_%s+_mapped.bam', CFG.tracks_dir, CFG.exp, CFG.genome_info.contig_names{c}), ...
                      sprintf('%s%s_%s+_el15_mm1_spliced.bam', CFG.tracks_dir, CFG.exp, CFG.genome_info.contig_names{c})};
end
CFG.tracks_max_intron_len = 1e9;
CFG.tracks_min_exon_len = -1;
CFG.tracks_max_mismatches = CFG.read_len;

%%%%% result directory %%%%%
date_exp = datestr(now,'yyyy-mm-dd');
%date_exp = datestr(now,'yyyy-mm-dd_HHhMM');
CFG.out_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/rquant/%s_%s/', CFG.exp, date_exp);
if ~exist(CFG.out_dir ,'dir'),
  [s m mid] = mkdir(CFG.out_dir);
  assert(s);
end

%%%%% optimizer %%%%%
CFG.optimizer = 'mosek';

%%%%% pre-learned profiles %%%%%
CFG.profiles_fn = '';

%%%%% rproc settings for rquant subjobs %%%%%
CFG.use_rproc = 1; % 1: cluster submission or 0: locally
if CFG.use_rproc,
  CFG.rproc_num_jobs              = 20;
  CFG.rproc_memreq                = 4000;
  CFG.rproc_par.priority          = 8;
  CFG.rproc_par.resubmit          = 3;
  CFG.rproc_par.mem_req_resubmit  = [8000 12000 20000];
  CFG.rproc_par.time_req_resubmit = [36*60 70*60 90*60];
  CFG.rproc_par.express           = 0;
  CFG.rproc_par.immediately_bg    = 0;
  CFG.rproc_par.immediately       = 0;
  CFG.rproc_par.arch              = 64;
  CFG.rproc_par.identifier        = '';
  CFG.rproc_par.verbosity         = 0;
  CFG.rproc_time                  = 5*60; % mins
end

CFG

%%%%% rproc settings for main job %%%%%
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
