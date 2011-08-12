addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% experiment %%%%%
CFG.organism = 'celegans';
CFG.exp = {'elegans.1_reads_0.fl.fastq.3a', 'elegans.2_reads_0.fl.fastq.3a', 'elegans.4_reads_0.fl.fastq.ffx.3a', 'elegans.5_reads_0.fl.fastq.ffx.3a'};
PAR.CFG.read_len = 75;

%%%%% tracks, repeats, genes, genome info %%%%% 
CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/celegans';
PAR.CFG.repeats_fn = '';
PAR.CFG.correct_intervals = 1;
PAR.anno_dir = sprintf('%s/annotation', CFG.base_dir);
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

for e = 1:length(CFG.exp),
  PAR.track = sprintf('%s/tracks/%s.bam', CFG.base_dir, CFG.exp{e});
  %%%%% result directory %%%%%
  date_exp = datestr(now,'yyyy-mm-dd');
  PAR.output_dir = sprintf('%s/rquant/%s_%s', CFG.base_dir, CFG.exp{e}, date_exp);
  PAR.output_file = sprintf('%s/%s_rquant.gff3', PAR.output_dir, CFG.exp{e});
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
    rproc_par.identifier = sprintf('rq.%s-', CFG.organism(1:2));
    fprintf(1, 'Submitting job %s to cluster\n', rproc_par.identifier);
    job = rproc('rquant_rproc', PAR, rproc_memreq, rproc_par, rproc_time);
    pause(10);
  else
    rquant_rproc(PAR);
  end
end