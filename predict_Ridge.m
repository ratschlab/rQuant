%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2007-2010 Georg Zeller, Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2007-2010 Max Planck Society
%

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


% generate numerical kmer vectors from sequences
if CFG.VERBOSE>1, tic; end
X_num = seq_2_kmers(X, CFG.RR.order);
if CFG.VERBOSE>1, fprintf(1, 'Feature generation for Ridge regression took %.1fs.\n', toc); end

N = size(X_num,2);
Y_pred = nan(1,N);
for n = 1:N,
  Y_pred(n) = w' * X_num(:,n);
end
assert(~any(isnan(Y_pred)));