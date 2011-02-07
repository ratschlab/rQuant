%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function [coverage excluded_reads ok intron_list read_starts] = get_coverage_per_read(CFG, gene, reverse_ret)
% [coverage excluded_reads ok intron_list read_starts] = get_coverage_per_read(CFG, gene, reverse_ret)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc.
% reverse_ret: reverse exon_count for minus strand
%
% -- output --
% coverage: matrix of exonic positions x reads
% excluded_reads: reads excluded per transcript
% ok: indicates success of file parsing
% intron_list: nx4 list of introns (intron start, intron stop, confirmation, strand)
% read_starts: vector of number of reads starting at the given exonic positions 


MIN_COV = 0.5;
ok = 1;
coverage = [];
excluded_reads = [];
intron_list = zeros(2, 0);
read_starts = zeros(gene.exonic_len, 1);

collapse = double(~CFG.paired & nargout<=4);

if nargin<3
  reverse_ret = 0;
end

if ~isfield(CFG, 'both_strands')
  CFG.both_strands = 0;
end

if CFG.both_strands
  strand = '0';
else
  strand = gene.strand;
end

if ~isfield(CFG, 'tracks_max_intron_len')
  CFG.tracks_max_intron_len = 1e9;
end

if ~isfield(CFG, 'tracks_min_exon_len')
  CFG.tracks_min_exon_len = -1;
end

if ~isfield(CFG, 'tracks_max_mismatches')
  CFG.tracks_max_mismatches = CFG.read_len;
end

win = CFG.read_len;
eidx = [max(gene.eidx(1)-win,1):gene.eidx(1)-1, gene.eidx, gene.eidx(end)+1:min(gene.eidx(end)+win,CFG.chr_len(gene.chr_num))];
win_size = length(max(gene.eidx(1)-win,1):gene.eidx(1)-1);

for f = 1:length(CFG.tracks_fn{gene.chr_num}),
  fname = CFG.tracks_fn{gene.chr_num}{f};
  if ~exist(fname, 'file'),
    warning('BAM file %s does not exist', fname);
  end
  try
    if nargout>3
      [coverage_idx_tmp{f}, intron_list_tmp] = get_reads(fname, gene.chr, eidx(1), eidx(end), strand, collapse, 1000, CFG.tracks_max_intron_len, CFG.tracks_min_exon_len, CFG.tracks_max_mismatches);
    else
      [coverage_idx_tmp{f}] = get_reads(fname, gene.chr, eidx(1), eidx(end), strand, collapse, 1000, CFG.tracks_max_intron_len, CFG.tracks_min_exon_len, CFG.tracks_max_mismatches);
    end
  catch
    warning('get_reads failed');
    intron_list = intron_list';
    ok = 0;
    return;
  end
  if exist('intron_list_tmp', 'var')
    if ~isempty(intron_list_tmp),
      intron_list = [intron_list intron_list_tmp];
    end
  end
end

% process coverage: convert to exonic position indices
coverage_idx = [coverage_idx_tmp{:}];
if ~collapse
  coverage_idx = unique(coverage_idx', 'rows')'; % no overlapping reads
  if ~isempty(coverage_idx)
    coverage = sparse(coverage_idx(1,:), coverage_idx(2,:), 1, max(coverage_idx(1,:)), eidx(end)-eidx(1)+1)';
  else
    coverage = sparse([], [], 1, eidx(end)-eidx(1)+1, 0);
  end
  coverage = coverage(eidx(win_size+1:win_size+gene.exonic_len)-eidx(1)+1, :);
  % no overlapping reads
  assert(~any(any(full(coverage>1))));
else
  coverage = sum(coverage_idx, 1);
  coverage = sparse(coverage(eidx(win_size+1:win_size+gene.exonic_len)-eidx(1)+1)');
end
  
% process intron list (1: intron start, 2: intron stop, 3: confirmation, 4: strand)
if nargout>3
  intron_list = [intron_list', zeros(size(intron_list,2), 1), (gene.strand=='-')*ones(size(intron_list,2), 1)];
  intron_list_unique = unique(intron_list, 'rows');
  for n = 1:size(intron_list_unique,1),
    intron_list_unique(n,3) = sum(intron_list_unique(n,1)==intron_list(:,1) & ...
                                  intron_list_unique(n,2)==intron_list(:,2) & ...
                                  intron_list_unique(n,4)==intron_list(:,4));
  end
  intron_list = intron_list_unique;
  clear intron_list_unique;
end

% process read starts: count reads starting at each exonic position
if nargout>4
  for c = 1:size(coverage, 2),
    fidx = find(coverage(:,c)~=0, 1, 'first');
    if ~isempty(fidx),
      if fidx==1 & sum(coverage(:,c),1)<CFG.read_len, continue; end
      read_starts(fidx) = read_starts(fidx) + 1;
    end
  end
end

% determine set of excluded reads
if CFG.paired
  offset = gene.start-1;
  eidx = 1:length(gene.eidx);
  sum_all = sum(coverage,1);
  for t = 1:length(gene.transcripts),
    tidx = [];
    for e = 1:size(gene.exons{t},1),
      tidx = [tidx, gene.exons{t}(e,1):gene.exons{t}(e,2)];
      tidx = unique(tidx);
    end
    tidx = unique(tidx)-offset;
    assert(all(tidx>=0));
    % transcript indices in relative exonic coordinates
    [tmp idx1 idx2] = intersect(tidx, gene.eidx-offset);
    assert(isequal(tidx, gene.eidx(idx2)-offset));
    tidx = idx2;
    % exclude reads that do not cover at least MIN_COV of all transcript positions
    excluded_reads{t} = find(sum(coverage(tidx,:),1)./sum_all<MIN_COV);
  end
else
 for t = 1:length(gene.transcripts),
   excluded_reads{t} = [];
 end
end

% collapse coverage if single reads are not required
if ~CFG.paired
  coverage = sum(coverage, 2);
end

% reverse for minus strand
if reverse_ret && gene.strand=='-'
  rev_idx = size(coverage,1):-1:1;
  coverage = coverage(rev_idx,:);
end
