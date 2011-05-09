function [genes num_del] = sanitise_genes(genes, CFG)
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
%   Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


% cell to struct
genes_cell = genes;
clear genes;
for g = 1:length(genes_cell), 
  gene = genes_cell{g};
  for e = 1:length(gene.exons)
    gene.exons{e} = double(gene.exons{e});
  end    
  gene.exons = reshape(gene.exons, 1, length(gene.exons));
  gene.id = double(gene.id);
  gene.start = double(gene.start);
  gene.stop = double(gene.stop);
  genes(g) = gene;
end
clear genes_cell;

% add exonic length
% initialise expression bins
% initialise transcript length bins
if ~isfield(genes, 'chr_num'),
  genes(1).chr_num = 0;
  chr = unique({genes.chr});
end
genes(1).expr_bin = [];
genes(1).transcript_len_bin = [];
genes(1).eidx = [];
genes(1).exonic_len = 0;
del_idx = false(1, length(genes));
for g = 1:length(genes),
  if length(genes(g).exons)==0,
    del_idx(g) = true;
    continue;
  end
  genes(g).chr_num = strmatch(genes(g).chr, chr, 'exact');
  eidx = [];
  min_start = inf; max_stop = -inf;
  % closed intervals
  if isequal(CFG.gene_source, 'annotation')
    genes(g).strands = repmat(genes(g).strand, 1, length(genes(g).transcripts));
  end
  if genes(g).start>genes(g).stop || genes(g).start<1 || genes(g).stop<1
    del_idx(g) = true;
  end
  genes(g).expr_bin = ones(1, length(genes(g).transcripts));
  genes(g).transcript_len_bin = ones(1, length(genes(g).transcripts));
  for t = 1:length(genes(g).transcripts),
    tidx = [];
    % closed intervals
    if isequal(CFG.gene_source, 'annotation')
      assert(isfield(genes, 'strands'));
      if genes(g).strands(t)=='+'
        genes(g).exons{t}(:,2) = genes(g).exons{t}(:,2)-1;
      else
        genes(g).exons{t}(:,1) = genes(g).exons{t}(:,1)+1;
      end
    end
    min_start = min(min(genes(g).exons{t}(:,1)), min_start);
    max_stop = max(max(genes(g).exons{t}(:,2)), max_stop);
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
  if isequal(CFG.gene_source, 'annotation')
    if (min_start~=genes(g).start)
      if ~(min_start==genes(g).start+1)
        %warning('annotated gene start does not agree with annotated transcript starts');
      end
    end
    if (max_stop~=genes(g).stop)
      if ~(max_stop==genes(g).stop-1)
        %warning('annotated gene stop does not agree with annotated transcript ends');
      end
    end
  end
  genes(g).start = min_start;
  genes(g).stop = max_stop;
  assert(genes(g).start<=genes(g).stop);
  genes(g).eidx = unique(eidx);
  genes(g).exonic_len = length(genes(g).eidx);
end
genes(del_idx) = [];
num_del = sum(del_idx);