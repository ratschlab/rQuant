addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% directories from which to load read data and genes %%%%%
CFG.organism = 'human'; % 'elegans' or 'human'
CFG.exp = 'TAlioto_full';
CFG.read_len = 75;
CFG.gene_source = 'annotation';

switch CFG.organism
 case 'human'
  CFG.tracks_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/eval/current_submissions/human_laneFiltered/';
  CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/human/hg19/repeat_masker/tracks/';
 case 'elegans'
  CFG.tracks_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/eval/current_submissions/elegans/';
  CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/repeat_masker/tracks/';
end

%%%%% genes %%%%%
switch CFG.organism
 case 'human'
  CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_nanostring/rgasp/nanostr_genes_HepG2.mat';
 case 'elegans'
  CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_nanostring/rgasp/nanostr_genes_C_elegans_L3.mat';
end
   
%%%%% genome info %%%%%
switch CFG.organism
 case 'human'
  CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/human/hg19/hg19.gio/genome.config');
 case 'elegans'
  CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/elegans/elegans.gio/genome.config');
end
% contig length
for c = 1:length(CFG.genome_info.flat_fnames),
  d = dir(CFG.genome_info.flat_fnames{c});
  CFG.chr_len(c) = d(1).bytes;
end

%%%%% alignments %%%%%
for c = 1:length(CFG.genome_info.flat_fnames),
  CFG.tracks_fn{c} = {sprintf('%s%s/sanitized.unique.region_filtered.bam', CFG.tracks_dir, CFG.exp)};
end
switch CFG.organism
 case 'human'
  CFG.tracks_max_intron_len = 200000; 
 case 'elegans'
  CFG.tracks_max_intron_len = 20000;
end
CFG.tracks_min_exon_len = 15;
CFG.tracks_max_mismatches = 0;
max_mismatches = [0 1 2];

%%%%% result directory %%%%
CFG.out_dir = '';

%%%%% optimizer %%%%%
CFG.optimizer = 'cplex';

%%%%% number of iterations (1: no profile learning) %%%%%
CFG.max_iter = 1;

%%%%% pre-learned profiles %%%%%
CFG.profiles_fn = '';

%%%%% rproc settings for rquant subjobs %%%%%
CFG.use_rproc = 0; % 1: cluster submission or 0: locally
if CFG.use_rproc,
  CFG.rproc_num_jobs              = 100;
  CFG.rproc_memreq                = 4000;
  CFG.rproc_par.priority          = 20;
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
rproc_memreq                = 7000;
rproc_par.priority          = 8;
rproc_par.express           = 0;
rproc_par.immediately_bg    = 0;
rproc_par.immediately       = 0;
rproc_par.arch              = 64;
rproc_par.identifier        = '';
rproc_par.verbosity         = 0;
rproc_time                  = 72*60; % mins
rproc_par.envstr            = 'export MOSEKLM_LICENSE_FILE=/fml/ag-raetsch/share/software/mosek/6/licenses/mosek.lic; export LD_LIBRARY_PATH=/fml/ag-raetsch/share/software/mosek/6/tools/platform/linux64x86/bin';


for m = 1:length(max_mismatches),
  CFG.tracks_max_mismatches = max_mismatches(m);
  CFG.out_fn = strrep(CFG.tracks_fn{1}{1}, '.bam', sprintf('.mm_%i.mat', CFG.tracks_max_mismatches));
  if 1
    rproc_par.identifier = sprintf('rq.%s-', CFG.organism(1:2));
    fprintf(1, 'Submitting job %i (%s) to cluster\n', m, rproc_par.identifier);
    job(m) = rproc('rquant', CFG, rproc_memreq, rproc_par, rproc_time);
  else
    rquant(CFG);
  end
end