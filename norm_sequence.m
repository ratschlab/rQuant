%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function norm_cov = norm_sequence(CFG, X)
% norm_cov = norm_sequence(CFG, X)
%
% -- input --
% CFG: configuration struct
% X: input data (k-mer dim x number of positions) 
%
% -- output --
% norm_cov: predicted normalised coverage

if ~any(isnan(CFG.RR.seq_norm_weights)),
  Y_pred = predict_Ridge(CFG, X, CFG.RR.seq_norm_weights);
  num_read_starts_pred = exp(Y_pred);
else
  num_read_starts_pred = ones(1, size(X,2));
end

norm_cov = zeros(length(num_read_starts_pred), 1);
for n = 1:length(num_read_starts_pred),
  idx = n:min(n+CFG.read_len-1, length(num_read_starts_pred));
  norm_cov(idx, 1) = norm_cov(idx) + num_read_starts_pred(n);
end
norm_cov = norm_cov / mean(norm_cov);