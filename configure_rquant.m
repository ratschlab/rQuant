%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

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

%%%%% transcript weight optimisation
% method to determine transcript weights 
CFG.method = 'pos'; % 'pos' or 'seg'
CFG.paired = 0;
% regularisation strength in transcript weight optimisation
CFG.C1 = 1;
CFG.C1_set = [0.001];
CFG.C1_loss_frac_target = 0.3;

%%%%% sequence bias normalisation
CFG.norm_seqbias = 0;
CFG.RR.seq_norm_weights = [];
CFG.RR.half_win_size = 20;
CFG.RR.num_train_frac = 0.8;
CFG.RR.order = 2;
CFG.RR.lambda = 1e-2;

%%%%% profile learning
% enables loading of profiles from CFG.profiles_fn
CFG.load_profiles = 1;
% number of plifs for profile functions
CFG.num_plifs = 50;
% maximal number of positions to be considered at both transcript ends
CFG.max_side_len = 5000;
% bins for different expression levels
exr = [-1 inf];
CFG.expr_ranges = round([exr(1:end-1)'+1 exr(2:end)']);
% bins for different transcript lengths % prctile(tlen,10) prctile(tlen,90)
tlr = [0 649 1008 1379 1977 inf];
%tlr = [1 750 921 1092 1263 1434 1605 1776 1947 2118 2289 2461 2632 2803 2974 3145 3316 3487 3658 3829 4000 Inf];
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
% will be weighted by the number of positions in profile_genes
CFG.C2.tau   = 100;
CFG.C2.kappa = 1;
CFG.C2.theta = 10;
% more output to stdout
CFG.VERBOSE = 1; % 0: no output, 1: more output, 2: debug output