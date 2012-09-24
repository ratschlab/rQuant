function X = gen_sequence_features(CFG, genes, tscp_idx)
% GEN_SEQUENCE_FEATURES   Generates sequence features for a set of genes.
%
%   X = gen_sequence_features(CFG, genes, tscp_idx)
%
%   -- input --
%   CFG:      configuration struct
%   genes:    struct defining genes with start, stops, exons etc.
%   tscp_idx: index of transcript for each gene
%
%   -- output --
%   X:        matrix of positional substring occurence (number of k-mers x positions)
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%

addpath(CFG.paths)

P = sum([genes.exonic_len]); % number of positions
S = 0;
for o = 1:CFG.seq.order,
  S = S + (CFG.seq.half_win_size*2-o+1) * 4^o;
end
X = sparse(S, P);
ci = 0;
if CFG.VERBOSE>1, fprintf(1, 'Generating features for sequence normalisation...\n'); tic; end
tmp_VERBOSE = CFG.VERBOSE;
CFG.VERBOSE = 0;
for g = 1:length(genes),
  if tmp_VERBOSE>1, fprintf(1, '%i\r', g); end
  exons = genes(g).exons{tscp_idx(g)};
  if ~isfield(genes(g), 'strands') || length(genes(g).strands)<tscp_idx(g),
    genes(g).strands(tscp_idx(g)) = genes(g).strand;
  end
  exons(1,1) = exons(1,1) - CFG.seq.half_win_size; 
  exons(end,2) = exons(end,2) + CFG.seq.half_win_size;
  seq = upper(load_genomic(genes(g).chr, genes(g).strands(tscp_idx(g)), exons(:,1), exons(:,2), CFG.genome_info, 0));
  % generate numerical kmer vectors from sequences
  tmp_X = seq_2_kmers(seq, CFG.seq.order, 2*CFG.seq.half_win_size);
  assert(size(tmp_X,1)==S);
  assert(size(tmp_X,2)<=genes(g).exonic_len);
  X(:, ci+[1:size(tmp_X,2)]) = tmp_X;
  clear tmp_X;
  ci = ci + genes(g).exonic_len;
end
CFG.VERBOSE = tmp_VERBOSE;
if CFG.VERBOSE>1, fprintf(1, 'Took %.1fs.\n', toc); end
