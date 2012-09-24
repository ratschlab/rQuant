function Y = gen_sequence_targets(CFG, genes, profile_weights)
% GEN_SEQUENCE_TARGETS   Generates target values for a set of genes.
%
%   Y = gen_sequence_targets(CFG, genes, profile_weights)
%
%   -- input --
%   CFG:             configuration struct
%   genes:           struct defining genes with start, stops, exons etc.
%   profile_weights: weights of profile functions
%
%   -- output --
%   Y:               vector of log ratio of observed and expected read starts
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%


P = sum([genes.exonic_len]); % number of positions
Y = sparse(1, P); % target values
ci = 0;
if CFG.VERBOSE>1, fprintf(1, 'Generating target values sequence normalisation...\n'); tic; end
tmp_VERBOSE = CFG.VERBOSE;
CFG.VERBOSE = 0;
for g = 1:length(genes),
  if tmp_VERBOSE>1, fprintf(1, '%i\r', g); end
  assert(length(genes(g).transcripts)==1);
  t = 1;
  assert(length(genes(g).num_read_starts)==sum(genes(g).exons{t}(:,2)-genes(g).exons{t}(:,1)+1));
  num_read_starts = genes(g).num_read_starts;
  % correct for profiles
  if isfield(genes, 'strands') && length(genes(g).strands)==length(genes(g).transcripts)
    strand_str = genes(g).strands(t);
  else
    strand_str = genes(g).strand;
  end
  if strand_str=='+'
    rev_idx = 1:size(profile_weights,1);
  else
    rev_idx = size(profile_weights,1):-1:1;
  end
  [feat feat_val feat_val_next] = gen_exon_features(CFG, genes(g), t);
  fidx = find(feat_val>0);
  feat = feat(fidx,:); feat_val = feat_val(fidx,:); feat_val_next = feat_val_next(fidx,:);
  exon_mask = gen_exon_mask(profile_weights(rev_idx,:), genes(g).transcript_len_bin(t), feat, feat_val, feat_val_next, [1:length(fidx)]', ones(length(fidx),1));
  num_read_starts = num_read_starts./exon_mask;
  % ratio of observed and expected read starts
  num_read_starts_norm = length(num_read_starts) * num_read_starts / sum(num_read_starts);
  tmp_Y = log(num_read_starts_norm + 1e-5)';
  Y(1, ci+[1:size(tmp_Y,2)]) = tmp_Y;
  clear tmp_Y;
  ci = ci + genes(g).exonic_len;
end
CFG.VERBOSE = tmp_VERBOSE;
if CFG.VERBOSE>1, fprintf(1, 'Took %.1fs.\n', toc); end