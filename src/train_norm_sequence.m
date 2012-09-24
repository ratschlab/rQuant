function [seq_weights, seq_feat, Q1, Q2] = train_norm_sequence(CFG, genes, seq_feat, profile_weights)
% TRAIN_NORM_SEQUENCE   Trains Ridge regression for sequence normalisation.
%
%   [seq_weights Q1 Q2] = train_norm_sequence(CFG, genes)
%
%   -- input --
%   CFG:         configuration struct
%   genes:       struct defining genes with start, stops, exons etc. 
%
%   -- output --
%   seq_weights: weights for sequence normalisation
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%


seq_weights = nan; 

if CFG.VERBOSE>0, fprintf(1, '\nStarting training for sequence normalisation...\n'); tic; end
%%%%% pre-processing %%%%%
if ~exist('seq_feat', 'var') || (exist('seq_feat', 'var') && isempty(seq_feat))
  tscp_idx = ones(1, length(genes));
  seq_feat = gen_sequence_features(CFG, genes, tscp_idx);
end
seq_target = gen_sequence_targets(CFG, genes, profile_weights);
assert(size(seq_feat,2)==size(seq_target,2));
fidx = find(~isnan(seq_target));
seq_feat = seq_feat(:,fidx);
seq_target = seq_target(:,fidx);
% split up data into training and test set
ridx = randperm(size(seq_feat,2));
num_train = floor(length(ridx)*CFG.seq.num_train_frac);
if num_train < 10^4
  if CFG.VERBOSE>0, fprintf(1, 'insufficient number of training examples (%i) for Ridge regression\n', num_train); end
  return;
end
X_train = seq_feat(:,ridx(1:num_train));
Y_train = seq_target(:,ridx(1:num_train));
X_eval = seq_feat(:,ridx(num_train+1:end));
Y_eval = seq_target(:,ridx(num_train+1:end));

%%%%% train weights of Ridge regression %%%%%
if CFG.VERBOSE>0, fprintf(1, '%i training examples for Ridge regression\n', size(X_train,2)); end
seq_weights = train_Ridge(X_train, Y_train, CFG.seq.lambda);

%%%%% performance measure %%%%%
if nargout>2
  % training set
  Y_pred_train = predict_Ridge(X_train, seq_weights);
  [Q1.TR Q2.TR] = mean_var_coeff(Y_pred_train, Y_train);
  % test set
  if CFG.VERBOSE>0, fprintf(1, '%i test examples for Ridge regression\n', size(X_eval,2)); end
  Y_pred_eval = predict_Ridge(X_eval, seq_weights);
  [Q1.TE Q2.TE] = mean_var_coeff(Y_pred_eval, Y_eval);
  if CFG.VERBOSE>0,
    fprintf(1, 'relative variability (Q1): training %.3f, test %.3f\n', Q1.TR, Q1.TE);
    fprintf(1, 'relative variability (Q2): training %.3f, test %.3f\n', Q2.TR, Q2.TE);
  end
end
if CFG.VERBOSE>0, fprintf(1, 'Took %.1fs.\n', toc); end