addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% directories from which to load read data and genes %%%%%
CFG.organism = 'elegans'; % 'arabidopsis' 'human'
CFG.exp = 'fs_strong_bias'; % 'fs_strong_bias_seq_bias'
CFG.gene_source = 'annotation';

%%%%% tracks, repeats, genes, genome info %%%%% 
switch CFG.organism
 case 'arabidopsis'
  %CFG.tracks_dir = '/fml/ag-raetsch/home/bohnert/tmp/ara_test/';
  %CFG.repeats_fn = '';
  %CFG.gene_fn = '/fml/ag-raetsch/home/bohnert/tmp/ara_test/genes_ara.mat';
  %CFG.genome_info = init_genome('/fml/ag-raetsch/home/bohnert/tmp/ara_test/genome_tair7.config');
  CFG.tracks_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10/tracks/';
  CFG.repeats_fn = '';
  CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10/run_2011-03-28/genes_expr.mat';
  CFG.genome_info = init_genome('/fml/ag-raetsch/share/databases/genomes/A_thaliana/arabidopsis_tair10/annotations/genome.gio/genome.config');
 case 'elegans'
  CFG.tracks_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/tracks/';
  CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/repeat_masker/tracks/';
  %CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/run_2011-01-26/genes_expr.mat';
  CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/run_2011-07-06/genes_expr.mat';
  CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/elegans/elegans.gio/genome.config');
 case 'human'
  CFG.tracks_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19/tracks/';
  CFG.repeats_fn = '';
  CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19/run_2011-05-02/genes_expr.mat';
  CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/human/hg19/hg19.gio/genome.config');
end
% contig length
for c = 1:length(CFG.genome_info.flat_fnames),
  d = dir(CFG.genome_info.flat_fnames{c});
  CFG.chr_len(c) = d(1).bytes;
end

%%%%% alignments %%%%%
CFG.samtools_dir = '/fml/ag-raetsch/share/software/samtools/';
for c = 1:length(CFG.genome_info.flat_fnames),
  CFG.tracks_fn{c} = {sprintf('%s%s.bam', CFG.tracks_dir, CFG.exp)};
end
CFG.tracks_max_intron_len = 1e9;
CFG.tracks_min_exon_len = -1;
CFG.tracks_max_mismatches = 1e3;

%%%%% output files %%%%%
CFG.write_gff = 0;
CFG.write_density_model = 0;

%%%%% enables profile learning %%%%%
CFG.learn_profiles = 1;

%%%%% pre-learned profiles %%%%%
CFG.profiles_fn = '';

%%%%% regularisation strengths %%%%%
C_I = 100;%[10^0 10^1 10^2];
C_F = 100;%[10^1 10^2 10^3 10^4];
C_N = 10;%[10^0 10^1 10^2];

%%%%% rproc settings for rquant subjobs %%%%%
CFG.use_rproc = 0; % 1: cluster submission or 0: locally
if CFG.use_rproc,
  CFG.rproc_num_jobs              = 50;
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

run_local = 1;

for s = 1:length(C_I),
  CFG.C_I = C_I(s);
  for f = 1:length(C_F),
    CFG.C_F = C_F(f);
    for n = 1:length(C_N),
      CFG.C_N = C_N(n);
      exp_str = sprintf('s%if%in%i', CFG.C_I, CFG.C_F, CFG.C_N);
      %%%%% result directory %%%%%
      %date_exp = datestr(now,'yyyy-mm-dd');
      date_exp = datestr(now,'yyyy-mm-dd_HHhMM');
      switch CFG.organism
       case 'arabidopsis'
        CFG.out_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10/rquant/%s_%s_%s/', CFG.exp, date_exp, exp_str);
       case 'elegans'
        CFG.out_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/rquant/%s_%s_%s/', CFG.exp, date_exp, exp_str);
       case 'human'
        CFG.out_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19/rquant/%s_%s_%s/', CFG.exp, date_exp, exp_str);
      end
      if ~exist(CFG.out_dir ,'dir'),
        [s m mid] = mkdir(CFG.out_dir);
        assert(s);
      end
      CFG.profiles_fn = sprintf('%s/profiles.mat', CFG.out_dir);
      %%%%% rproc settings for main job %%%%%
      if ~run_local
        rproc_memreq                = 2000;
        rproc_par.priority          = 8;
        rproc_par.express           = 0;
        rproc_par.immediately_bg    = 0;
        rproc_par.immediately       = 0;
        rproc_par.arch              = 64;
        rproc_par.identifier        = '';
        rproc_par.verbosity         = 0;
        rproc_time                  = 72*60;
        rproc_par.identifier = sprintf('rq.%s%s-', CFG.organism(1:2), exp_str);
        fprintf(1, 'Submitting job %s to cluster\n', rproc_par.identifier);
        job = rproc('rquant', CFG, rproc_memreq, rproc_par, rproc_time);
      else
        rquant(CFG);
      end
    end
  end
end