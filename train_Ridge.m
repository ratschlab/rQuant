function [w TR] = train_Ridge(CFG, X, Y)
% [w TR] = train_Ridge(CFG, X, Y)
%
% -- input --
% CFG: configuration struct
% X: input data (sequences of length 2*CFG.RR.half_win_size)
% Y: target values (expression values)
%
% -- output --
% w: weight vector of trained Ridge regression
% TR: training errors of trained Ridge regression


% generate numerical kmer vectors from sequences
if CFG.VERBOSE>0, tic; end
X_num = seq_2_kmers(X, CFG.RR.order);
if CFG.VERBOSE>0, fprintf(1, 'Feature generation for Ridge regression took %.1fs.\n', toc); end

if CFG.VERBOSE>0, tic; end
% regularizer
x_dim = size(X_num, 1);
R = CFG.RR.lambda * eye(x_dim, x_dim);
% number of training examples
N = size(X_num,2);
Q = X_num * X_num';
b = X_num * Y';
w = (R + Q)\b;
assert(~any(isnan(w)));
if CFG.VERBOSE>0, fprintf(1, 'Training for Ridge regression took %.1fs.\n', toc); end

% training error
Y_pred = nan(1,N);
for n = 1:N,
  Y_pred(n) = w' * X_num(:,n);
end

if CFG.VERBOSE>0,
  % absolute variablity on training set
  TR.Q1 = mean(abs(Y_pred - Y)) / mean(abs(Y - median(Y)));
  % squared variablity on training set
  TR.Q2 = mean((Y_pred - Y).^2) / mean((Y - mean(Y)).^2);
end
