function CFG = configure_rquant(CFG)
% configure_rquant(CFG)
%
% -- input --
% CFG.organism
% CFG.exp
% CFG.read_len
% CFG.gene_source
% CFG.read_maps_select
% CFG.fname
% CFG.paths


%%%%% directories from which to load read data and genes %%%%%
CFG.base_dir = '../examples/';
CFG.read_maps_dir = sprintf('%stracks/', CFG.base_dir);
CFG.repeat_maps_dir = sprintf('%sannotations/%s/repeat_masker/tracks/', CFG.base_dir, CFG.organism);
CFG.repeat_maps_fn = CFG.repeat_maps_dir;
CFG.profiles_fn = '';
%CFG.profiles_fn = sprintf('%srquant/%s/%s/profiles.mat', CFG.base_dir, CFG.organism, CFG.exp);
CFG.samtools_dir = '/fml/ag-raetsch/share/software/samtools-0.1.7a/';

switch CFG.gene_source
 case 'annotation'
  CFG.gene_dir = sprintf('%sannotations/', CFG.base_dir);
  switch CFG.organism,
   case 'drosophila'
    CFG.gene_fn = '';
   case 'elegans'
    CFG.gene_fn = sprintf('%sgenes.mat', CFG.gene_dir);
   case 'human'
    CFG.gene_fn = '';
   otherwise
    error('unknown organism %s', organism);
  end
 case 'mtim'
  CFG.gene_dir = '';
  CFG.gene_fn = '';
 case 'mgene'
  CFG.gene_dir = '';
  CFG.gene_fn = '';  
 otherwise
  error('unknown gene source %s', CFG.gene_source);
end


%%%%% genome config %%%%%
switch CFG.organism,
 case 'drosophila'
  CFG.genome_info = '';
 case 'elegans'
  CFG.genome_info = init_genome(sprintf('%sgenomes/%s/genome.config', CFG.base_dir, CFG.organism));
 case 'human'
  CFG.genome_info = '';
 otherwise
  error('unknown organism %s', organism);
end
% contig length
for c = 1:length(CFG.genome_info.flat_fnames),
  d = dir(CFG.genome_info.flat_fnames{c});
  CFG.chr_len(c) = d(1).bytes;
end

for c = 1:length(CFG.genome_info.flat_fnames),
  CFG.introns_fn{c} = {sprintf('%s%s/%s/%s_%s+%s.introns', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select),
                       sprintf('%s%s/%s/%s_%s-%s.introns', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select)};
  CFG.read_maps_fn{c} = {sprintf('%s%s/%s/%s_%s+%s.bam', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select)};
  %CFG.read_maps_fn{c} = {sprintf('%s%s/%s/%s_%s+%s_spliced.bam', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select)};
  %CFG.read_maps_fn{c} = {sprintf('%s%s/%s/%s_%s+%s_mapped.bam', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select), ...
   %                      sprintf('%s%s/%s/%s_%s+%s_spliced.bam', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select), ...
    %                     sprintf('%s%s/%s/%s_%s-%s_spliced.bam', CFG.read_maps_dir, CFG.organism, CFG.exp, CFG.exp, CFG.genome_info.contig_names{c}, CFG.read_maps_select),};
end

%%%%% directory where to store results %%%%%
if isequal(CFG.gene_source, 'annotation')
  CFG.out_dir = sprintf('%srquant/%s/%s/rquant_result_%s/', CFG.base_dir, CFG.organism, CFG.exp, datestr(now,'yyyy-mm-dd'));
  %CFG.out_dir = sprintf('%srquant/%s/%s/rquant_result_%s/', CFG.base_dir, CFG.organism, CFG.exp, datestr(now,'yyyy-mm-dd_HHhMM'));
  if ~exist(CFG.out_dir ,'dir'),
    [s m mid] = mkdir(CFG.out_dir);
    assert(s);
  end
end


%%%%% rproc settings %%%%%
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


%%%%% rquant parameters %%%%%

% enables taking data for both strands together
if isequal(CFG.gene_source, 'annotation')
  CFG.both_strands = 1;
  %CFG.both_strands = 0;
else
  CFG.both_strands = 1;
end
% number of iterations (1: no profile learning)
if isequal(CFG.gene_source, 'annotation')
  %CFG.max_iter = 1;
  CFG.max_iter = 6;
else
  CFG.max_iter = 1;
end
% optimizer
CFG.optimizer = 'mosek';

%%%%% transcript weight optimisation
% method to determine transcript weights 
CFG.method = 'pos'; % 'pos' or 'seg'
CFG.paired = 0;
% regularisation strength in transcript weight optimisation
CFG.C1 = 1;
CFG.C1_set = [0.001];
CFG.C1_loss_frac_target = 0.3;

%%%%% sequence bias normalisation
CFG.norm_seqbias = 1;
CFG.RR.seq_norm_weights = [];
CFG.RR.half_win_size = 20;
CFG.RR.num_train_frac = 0.8;
CFG.RR.order = 2;
CFG.RR.lambda = 1e-2;

%%%%% profile learning
% enables loading of profiles from CFG.profiles_fn
CFG.load_profiles = 0;
% number of plifs for profile functions
CFG.num_plifs = 50;
% maximal number of positions to be considered at both transcript ends
CFG.max_side_len = 5000;
% bins for different expression levels
exr = [-1 inf];
CFG.expr_ranges = round([exr(1:end-1)'+1 exr(2:end)']);
% bins for different transcript lengths % prctile(tlen,10) prctile(tlen,90)
switch CFG.organism
 case 'drosophila'
  tlr = [0 1058 1645 2337 3382 inf];
 case 'elegans'
  tlr = [0 649 1008 1379 1977 inf];
 case 'human'
  tlr = [0 297 486 586 745 959 1451 2013 2738 4036 inf];
 otherwise
  error('unknown organism %s', organism);
end
CFG.transcript_len_ranges = round([tlr(1:end-1)'+1 tlr(2:end)']);
% bins for distances to closest intron
CFG.num_intron_plifs = 5;
% enables subsampling of data for learning profiles
CFG.subsample = 1;
% maximal number of examples for learning profiles
CFG.max_num_train_exm = 4e6 * 5;
% fraction of genes to be subsampled for learning profiles
CFG.subsample_frac_global = 1;
% fraction of profile_genes to be subsampled for learning profiles
CFG.subsample_frac = 0.10;
% regularisation strength in profile optimisation
CFG.C2.tau   = 100;
CFG.C2.kappa = 1;
CFG.C2.theta = 10;
CFG.C2.tau   = CFG.C2.tau*CFG.max_num_train_exm;
CFG.C2.kappa = CFG.C2.kappa*CFG.max_num_train_exm;
CFG.C2.theta = CFG.C2.theta*CFG.max_num_train_exm;
% more output to stdout
CFG.VERBOSE = 1; % 0: no output, 1: more output, 2: debug output


%%%%% paths dependent on configuration %%%%%
if CFG.use_rproc
  p = '~/svn/tools/rproc';
  CFG.paths = sprintf('%s:%s', p, CFG.paths);
  addpath(p);
end

p = '/fml/ag-raetsch/share/software/matlab_tools/cplex9';
CFG.paths = sprintf('%s:%s', p, CFG.paths);
addpath(p);
switch CFG.optimizer
 case 'cplex'
  p = '/fml/ag-raetsch/share/software/matlab_tools/cplex9';
  CFG.paths = sprintf('%s:%s', p, CFG.paths);
  addpath(p);
 case 'mosek'
  p = '/fml/ag-raetsch/share/software/mosek/5/toolbox/r2007a';
  CFG.paths = sprintf('%s:%s', p, CFG.paths);
  addpath(p);
  if CFG.use_rproc
    CFG.rproc_par.envstr = 'export MOSEKLM_LICENSE_FILE=/fml/ag-raetsch/home/bohnert/tmp/mosek.lic; export LD_LIBRARY_PATH=/fml/ag-raetsch/share/software/mosek/5/tools/platform/linux64x86/bin';
  end
 otherwise
  error('unknown optimizer %s', CFG.optimizer);
end
