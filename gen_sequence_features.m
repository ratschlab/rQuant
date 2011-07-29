function X = gen_sequence_features(CFG, gene, tscp_idx)
% GEN_SEQUENCE_FEATURES   Generates sequence features for a transcript.
%
%   X = gen_sequence_features(CFG, gene, tscp_idx)
%
%   -- input --
%   CFG:      configuration struct
%   gene:     struct defining a gene with start, stops, exons etc.
%   tscp_idx: index of transcript
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


exons = gene.exons{tscp_idx};
eidx = [];
for e = 1:size(exons,1),
  eidx = [eidx, exons(e,1):exons(e,2)];
end
eidx = unique(eidx);
[tmp idx1 idx2] = intersect(eidx, gene.eidx);

if ~isfield(gene, 'strands') || length(gene.strands)<tscp_idx,
  gene.strands(tscp_idx) = gene.strand;
end
exons(1,1) = exons(1,1) - CFG.seq.half_win_size; 
exons(end,2) = exons(end,2) + CFG.seq.half_win_size;
seq = upper(load_genomic(gene.chr, gene.strands(tscp_idx), exons(:,1), exons(:,2), CFG.genome_info, 0));

% generate numerical kmer vectors from sequences
X = seq_2_kmers(seq, CFG.seq.order, 2*CFG.seq.half_win_size);

