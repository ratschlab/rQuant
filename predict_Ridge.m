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
% X: input data (k-mer dim x number of positions)
% w: weight vector of trained Ridge regression
%
% -- output --
% Y_pred: predicted target values


N = size(X, 2);
Y_pred = nan(1, N);
for n = 1:N,
  Y_pred(n) = w' * X(:,n);
end
assert(~any(isnan(Y_pred)));