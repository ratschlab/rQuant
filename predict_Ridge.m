function Y_pred = predict_Ridge(CFG, X, w)
% Y_pred = predict_Ridge(CFG, X, w)
%
% -- input --
% CFG: configuration struct
% X: input data (sequences of length 2*CFG.RR.half_win_size)
% w: weight vector of trained Ridge regression
%
% -- output --
% Y_pred: predicted target values


if CFG.VERBOSE>0, fprintf(1, '%i test examples for Ridge regression\n', size(X,2)); end

% generate numerical kmer vectors from sequences
if CFG.VERBOSE>0, tic; end
X_num = seq_2_kmers(X, CFG.RR.order);
if CFG.VERBOSE>0, fprintf(1, 'Feature generation for Ridge regression took %.1fs.\n', toc); end

N = size(X_num,2);
Y_pred = nan(1,N);
for n = 1:N,
  Y_pred(n) = w' * X_num(:,n);
end
assert(~any(isnan(Y_pred)));
