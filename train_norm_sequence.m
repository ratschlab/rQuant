%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function w = train_norm_sequence(CFG, X, Y)
% w = train_norm_sequence(CFG, X, Y)
%
% -- input --
% CFG: configuration struct
% X: matrix of positional substring occurence (number of k-mers x positions)
% Y: vector of target values for regression
%
% -- output --
% w: weight vector of trained Ridge regression


w = nan;

ridx = randperm(size(X,2));
num_train = floor(length(ridx)*CFG.RR.num_train_frac);

if num_train < 10^4
  if CFG.VERBOSE>0, fprintf(1, 'insufficient number of training examples (%i) for Ridge regression\n', num_train); end
  %return;
end

X_train = X(:,ridx(1:num_train));
Y_train = Y(:,ridx(1:num_train));
X_eval = X(:,ridx(num_train+1:end));
Y_eval = Y(:,ridx(num_train+1:end));

% train weights of Ridge regression
if CFG.VERBOSE>0, fprintf(1, '%i training examples for Ridge regression\n', size(X_train,2)); end
[w TR] = train_Ridge(CFG, X_train, Y_train);
% test on omitted examples
if CFG.VERBOSE>0,
  fprintf(1, '%i test examples for Ridge regression\n', size(X_eval,2));
  Y_pred = predict_Ridge(CFG, X_eval, w);
  % absolute variablity on test set
  TE.Q1 = mean(abs(Y_pred - Y_eval)) / mean(abs(Y_eval - median(Y_eval)));
  % squared variablity on test set
  TE.Q2 = mean((Y_pred - Y_eval).^2) / mean((Y_eval - mean(Y_eval)).^2);
  fprintf(1, 'relative variability (Q1): training %.3f, test %.3f\n', TR.Q1, TE.Q1);
  fprintf(1, 'relative variability (Q2): training %.3f, test %.3f\n', TR.Q2, TE.Q2);
end