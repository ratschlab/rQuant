function coverage_norm = norm_sequence(CFG, gene, tscp_idx, seq_weights)
% NORM_SEQUENCE   Predicts normalised sequence profile from sequence.
%
%   coverage_norm = norm_sequence(CFG, gene, tscp_idx, seq_weights)
%
%   -- input --
%   CFG:           configuration struct
%   gene:          struct defining a gene with start, stops, exons etc.
%   tscp_idx:      index of transcript
%   seq_weights:   weights for sequence normalisation
%
%   -- output --
%   coverage_norm: normalised sequence profile 
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


% extend by read length at boundaries
gene.exons{tscp_idx}(1,1) = gene.exons{tscp_idx}(1,1)-CFG.read_len+1;
gene.exons{tscp_idx}(end,2) = gene.exons{tscp_idx}(end,2)+CFG.read_len-1;
gene.exonic_len = sum(gene.exons{tscp_idx}(:,2)-gene.exons{tscp_idx}(:,1)+1);
% generate features
X = gen_sequence_features(CFG, gene, tscp_idx);
assert(size(X,1)==size(seq_weights,1));
assert(size(X,2)==gene.exonic_len);
% predict target values
Y_pred = predict_Ridge(X, seq_weights);
num_read_starts_pred = exp(Y_pred);
% compute normalised sequence profile
coverage = zeros(1, length(num_read_starts_pred));
for n = 1:length(num_read_starts_pred),
  idx = n:min(n+CFG.read_len-1, length(num_read_starts_pred));
  coverage(idx) = coverage(idx) + num_read_starts_pred(n);
end
% cut back to actual boundaries
coverage = coverage(CFG.read_len:end-CFG.read_len+1);
coverage_norm = coverage./mean(coverage);

