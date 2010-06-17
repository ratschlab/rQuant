function CFG = configure_rquant(CFG)
% configure_rquant(CFG)
%
% -- input --
% CFG: configure struct with paths, data directories etc.
%
% --output --
% CFG: configure struct augmented by rQuant parameters


%%%%% rquant parameters %%%%%

% enables taking data for both strands together
if isequal(CFG.gene_source, 'annotation')
  CFG.both_strands = 1;
  %CFG.both_strands = 0;
else
  CFG.both_strands = 1;
end
% number of iterations (1: no profile learning)
CFG.max_iter = 6;

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
tlr = [0 649 1008 1379 1977 inf];
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