%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function [feat, tlen] = gen_exon_features(gene, t, num_bins, max_side_len, reverse_ret)
% [feat, tlen] = gen_exon_features(gene, t, num_bins)
%
% -- input --
% gene: struct defining a gene with start, stops, exons etc.
% t: index of transcript
% num_bins: number of bins for PLiFs
% max_side_len: maximal number of positions to be considered at both transcript parts
%
% -- output --
% feat: P x num_bins matrix of features for P exonic positions
% tlen: length of transcript

if nargin<5
  reverse_ret = 0;
end

offset = gene.start-1;
exons = gene.exons{t};

% transcript indices in relative gene coordinates for all exons

eidx = gene.eidx;
eidx = unique(eidx-offset);

% transcript indices in relative gene coordinates
tidx = [];
for e = 1:size(exons,1),
  tidx = [tidx, exons(e,1):exons(e,2)];
  tidx = unique(tidx);
end
tidx = unique(tidx)-offset;
assert(all(tidx>=0));
tlen = length(tidx);

% transcript indices in relative exonic coordinates
[tmp idx1 idx2] = intersect(tidx, eidx);
assert(isequal(tidx, eidx(idx2)));
tidx = idx2;

feat = zeros(length(eidx), num_bins);
lmt = linspace(0, sqrt(max_side_len), (num_bins/2)+1).^2;
lmt(end) = inf;

% left transcript part
dist = 1:ceil(length(tidx)*0.5);
for b = 1:(num_bins/2),
  fidx = find(lmt(b)<=dist & lmt(b+1)>dist); 
  feat(tidx(fidx), b) = 1;
end

% right transcript part
dist = 1:(length(tidx)-length(tidx)*0.5);
for b = 1:(num_bins/2),
  fidx = find(lmt(b)<=dist & lmt(b+1)>dist); 
  feat(tidx(length(tidx)-fidx+1), num_bins-b+1) = 1;
end

% reverse for minus strand
if reverse_ret && gene.strand=='-'
  rev_idx = size(feat,1):-1:1;
  feat = feat(rev_idx,:);
end

assert(all(sum(feat,2)==1|sum(feat,2)==0));
assert(size(feat,1)==length(eidx));
