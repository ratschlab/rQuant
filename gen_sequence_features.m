function [X Y Y_bg] = gen_sequence_features(CFG, gene, tscp_idx)
% [X Y] = gen_sequence_features(CFG, gene, tscp_idx)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc.
% tscp_idx: index of transcript
%
% -- output --
% X: matrix of positional substring occurence (number of k-mers x positions)
% Y: vector of target values for regression
% Y_bg: vector of averaged number of read starts (background)

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
exons(1,1) = exons(1,1) - CFG.RR.half_win_size; 
exons(end,2) = exons(end,2) + CFG.RR.half_win_size;
seq = upper(load_genomic(gene.chr, gene.strands(tscp_idx), exons(:,1), exons(:,2), CFG.genome_info, 0));

% generate numerical kmer vectors from sequences
X = seq_2_kmers(seq, CFG.RR.order, 2*CFG.RR.half_win_size);

% target values for regression
if nargout>1
  num_read_starts = gene.num_read_starts(idx2);
  assert(length(seq)==length(num_read_starts)+2*CFG.RR.half_win_size);
  %num_read_starts_norm = length(num_read_starts) * num_read_starts / sum(num_read_starts);
  Y_bg = repmat(log(median(num_read_starts) + 1e-10), 1, length(num_read_starts));
  Y = log(num_read_starts + 1e-10)'-Y_bg;
  assert(size(X,2)==length(Y));
end
