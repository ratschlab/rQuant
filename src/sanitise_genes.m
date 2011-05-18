function [genes parent_genes num_del num_merged] = sanitise_genes(genes, CFG)
% SANITISE_GENES   Adapts genometool gene struct to rQuant internal gene struct
%
%   [genes num_del] = sanitise_genes(genes, CFG)
%
%   -- input --
%   genes:   struct defining genes with start, stops, exons etc.
%   CFG:     configuration struct
%
%   -- output --
%   genes:   augmented by eidx, exonic_len
%   num_del: number of deleted genes
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


assert(isstruct(genes));
% add chr_num field
if ~isfield(genes, 'chr_num'),
  [genes.chr_num] = deal(0);
  chr = unique({genes.chr});
  for g = 1:length(genes),
    genes(g).chr_num = strmatch(genes(g).chr, chr, 'exact');
  end
end

parent_genes = genes;
clear genes;
% merge transcripts from overlapping loci
genes = merge_transcripts_by_colocation(parent_genes, {'transcripts', 'exons'}, 0, CFG.VERBOSE>1);
num_merged = length(parent_genes) - length(genes);

% add exonic length
% initialise expression bins
% initialise transcript length bins
[genes.expr_bin] = deal([]);
[genes.transcript_len_bin] = deal([]);
[genes.eidx] = deal([]);
[genes.exonic_len] = deal(0);
del_idx = false(1, length(genes));
add_strands = ~isfield(genes, 'strands');
for g = 1:length(genes),
  if length(genes(g).exons)==0,
    del_idx(g) = true;
    continue;
  end
  eidx = [];
  if add_strands
    genes(g).strands = repmat(genes(g).strand, 1, length(genes(g).transcripts));
  end
  if genes(g).start>genes(g).stop || genes(g).start<1 || genes(g).stop<1
    del_idx(g) = true;
  end
  genes(g).expr_bin = ones(1, length(genes(g).transcripts));
  genes(g).transcript_len_bin = ones(1, length(genes(g).transcripts));
  for t = 1:length(genes(g).transcripts),
    tidx = [];
    if any(genes(g).exons{t}(:,1)>genes(g).exons{t}(:,2)) || any(genes(g).exons{t}(:,1)<1) || any(genes(g).exons{t}(:,2)<1) || ...
      length(genes(g).exons{t}(:,1))~=length(unique(genes(g).exons{t}(:,1))) || length(genes(g).exons{t}(:,2))~=length(unique(genes(g).exons{t}(:,2))),
      del_idx(g) = true;
    end
    for e = 1:size(genes(g).exons{t},1),
      tidx = [tidx, genes(g).exons{t}(e,1):genes(g).exons{t}(e,2)];
      tidx = unique(tidx);
    end
    genes(g).transcript_length(t) = length(unique(tidx));
    if genes(g).transcript_length(t)>100000,
      del_idx(g) = true;
    end
    eidx = unique([eidx tidx]);
  end
  assert(genes(g).start<=genes(g).stop);
  genes(g).eidx = unique(eidx);
  genes(g).exonic_len = length(genes(g).eidx);
end
genes(del_idx) = [];
num_del = sum(del_idx);