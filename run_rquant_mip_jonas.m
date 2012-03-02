addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% experiment %%%%%
CFG.organism = 'arabidopsis';
CFG.exp = {{'genes_expr_weak_bias1.bam'}, {'genes_expr_weak_bias2.bam'}, ...
           {'genes_expr_weak_bias3.bam'}, {'genes_expr_weak_bias4.bam'}};
PAR.CFG.read_len = 75;

%%%%% tracks, repeats, genes, genome info %%%%% 
CFG.base_dir = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt10_rquant';
CFG.out_dir = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt10_rquant';
PAR.CFG.repeats_fn = '';
PAR.CFG.correct_intervals = 0;
PAR.anno_dir = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt10_rquant';
PAR.track = '';

%%%%% output files %%%%%
PAR.output_dir = '';
PAR.output_file = '';
PAR.CFG.write_gff = 0;
PAR.CFG.write_density_model = 0;
PAR.profiles_fn_out = '';

%%%%% profile learning %%%%%
% enables profile learning
PAR.learn_profiles = 1; % 0: no learning, 1: empirically estimated, 2: optimised
% pre-learned profiles
PAR.load_profiles = 0;
PAR.profiles_fn = '';
% regularisation strengths
PAR.CFG.C_I = 100;
PAR.CFG.C_F = 100;
PAR.CFG.C_N = 10;
PAR.CFG.C_PE = 3*10^2;

%%%%% paired-end %%%%%
PAR.CFG.paired = 0;

%%%%% sequence bias normalisation %%%%%
PAR.CFG.norm_seqbias = 0;

%%%%% rproc settings for rquant subjobs %%%%%
PAR.CFG.use_rproc = 1; % 1: cluster submission or 0: locally
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

for e = 1:length(CFG.exp),
  for f = 1:length(CFG.exp{e}),
    PAR.track{f} = sprintf('%s/%s.bam', CFG.base_dir, CFG.exp{e}{f});
  end
  %%%%% result directory %%%%%
  %date_exp = datestr(now,'yyyy-mm-dd');
  date_exp = datestr(now,'yyyy-mm-dd_HHhMM')
  PAR.output_dir = sprintf('%s/rquant/%s_%s', CFG.out_dir, CFG.exp{e}{1}, date_exp);
  PAR.output_file = sprintf('%s/%s_rquant.gff3', PAR.output_dir, CFG.exp{e}{1});
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
    rproc_memreq                = 5000;
    rproc_par.priority          = 8;
    rproc_par.express           = 0;
    rproc_par.immediately_bg    = 0;
    rproc_par.immediately       = 0;
    rproc_par.arch              = 64;
    rproc_par.identifier        = '';
    rproc_par.verbosity         = 0;
    rproc_time                  = 72*60;
    rproc_par.identifier = sprintf('rq.%s-%i', CFG.organism(1:2), e);
    fprintf(1, 'Submitting job %s to cluster\n', rproc_par.identifier);
    job = rproc('rquant_rproc', PAR, rproc_memreq, rproc_par, rproc_time);
    pause(60);
  else
    rquant_rproc(PAR);
  end
end
