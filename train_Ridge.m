%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2007-2010 Georg Zeller, Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2007-2010 Max Planck Society
%

function [w TR] = train_Ridge(CFG, X, Y)
% [w TR] = train_Ridge(CFG, X, Y)
%
% -- input --
% CFG: configuration struct
% X: input data (k-mer dim x number of positions)
% Y: target values (expression values)
%
% -- output --
% w: weight vector of trained Ridge regression
% TR: training errors of trained Ridge regression


if CFG.VERBOSE>0, tic; end
% regularizer
x_dim = size(X, 1);
R = CFG.RR.lambda * eye(x_dim, x_dim);
% number of training examples
N = size(X,2);
Q = X * X';
b = X * Y';
w = (R + Q)\b;
assert(~any(isnan(w)));
if CFG.VERBOSE>0, fprintf(1, 'Training for Ridge regression took %.1f s.\n', toc); end

% training error
Y_pred = nan(1,N);
for n = 1:N,
  Y_pred(n) = w' * X(:,n);
end

if CFG.VERBOSE>0,
  % absolute variablity on training set
  TR.Q1 = mean(abs(Y_pred - Y)) / mean(abs(Y - median(Y)));
  % squared variablity on training set
  TR.Q2 = mean((Y_pred - Y).^2) / mean((Y - mean(Y)).^2);
end
