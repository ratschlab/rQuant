addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% experiment %%%%%
CFG.organism = 'elegans'; % 'arabidopsis' 'human'
CFG.exp = 'fs_strong_bias'; % 'fs_strong_bias_seq_bias'
PAR.CFG.read_len = 75;

%%%%% tracks, repeats, genes, genome info %%%%% 
switch CFG.organism
 case 'arabidopsis'
  CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10';
  PAR.CFG.repeats_fn = '';
  PAR.CFG.correct_intervals = 1;
 case 'elegans'
  CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200';
  PAR.CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/repeat_masker/tracks';
  PAR.CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/elegans/elegans.gio/genome.config');
  PAR.CFG.correct_intervals = 1;
 case 'human'
  CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19';
  PAR.CFG.repeats_fn = '';
  PAR.CFG.correct_intervals = 1;
end
PAR.anno_dir = sprintf('%s/annotation', CFG.base_dir);
PAR.track = sprintf('%s/tracks/%s.bam', CFG.base_dir, CFG.exp);

%%%%% output files %%%%%
PAR.output_dir = '';
PAR.output_file = '';
PAR.CFG.write_gff = 0;
PAR.CFG.write_density_model = 0;
PAR.profiles_fn_out = '';

%%%%% profile learning %%%%%
% enables profile learning
PAR.learn_profiles = 2; % 0: no learning, 1: empirically estimated, 2: optimised
% pre-learned profiles
PAR.load_profiles = 0;
PAR.profiles_fn = '';
% regularisation strengths
C_I = 100; %[10^0 10^1 10^2 10^3];
C_F = [10^1 10^2 10^3 10^4];
C_N = [10^0 10^1 10^2];

%%%%% sequence bias normalisation %%%%%
PAR.CFG.norm_seqbias = 0;

%%%%% rproc settings for rquant subjobs %%%%%
PAR.CFG.use_rproc = 0; % 1: cluster submission or 0: locally
if PAR.CFG.use_rproc,
  PAR.CFG.rproc_num_jobs              = 50;
  PAR.CFG.rproc_memreq                = 4000;
  PAR.CFG.rproc_par.priority          = 8;
  PAR.CFG.rproc_par.resubmit          = 3;
  PAR.CFG.rproc_par.mem_req_resubmit  = [8000 12000 20000];
  PAR.CFG.rproc_par.time_req_resubmit = [36*60 70*60 90*60];
  PAR.CFG.rproc_par.express           = 0;
  PAR.CFG.rproc_par.immediately_bg    = 0;
  PAR.CFG.rproc_par.immediately       = 0;
  PAR.CFG.rproc_par.arch              = 64;
  PAR.CFG.rproc_par.identifier        = '';
  PAR.CFG.rproc_par.verbosity         = 0;
  PAR.CFG.rproc_time                  = 5*60;
end


run_local = 0;

for s = 1:length(C_I),
  PAR.CFG.C_I = C_I(s);
  for f = 1:length(C_F),
    PAR.CFG.C_F = C_F(f);
    for n = 1:length(C_N),
      PAR.CFG.C_N = C_N(n);
      %%%%% result directory %%%%%
      %date_exp = datestr(now,'yyyy-mm-dd');
      date_exp = datestr(now,'yyyy-mm-dd_HHhMM')
      switch CFG.organism
       case 'arabidopsis'
        PAR.output_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10/rquant/%s_%s', CFG.exp, date_exp);
       case 'elegans'
        PAR.output_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/rquant/%s_%s', CFG.exp, date_exp);
       case 'human'
        PAR.output_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19/rquant/%s_%s', CFG.exp, date_exp);
      end
      PAR.output_file = sprintf('%s/%s_rquant.gff3', PAR.output_dir, CFG.exp);
      if ~exist(PAR.output_dir ,'dir'),
        [s m mid] = mkdir(PAR.output_dir);
        assert(s);
      end
      PAR.profiles_fn = sprintf('%s/profiles.mat', PAR.output_dir);
      % save parameters
      fname = strrep(PAR.output_file, '.gff3', '.par');
      write_parameters(PAR, fname);
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
        rproc_par.identifier = sprintf('rq.%s-', CFG.organism(1:2));
        fprintf(1, 'Submitting job %s to cluster\n', rproc_par.identifier);
        job = rproc('rquant_rproc', PAR, rproc_memreq, rproc_par, rproc_time);
        pause(60);
      else
        rquant_rproc(PAR);
      end
    end
  end
end