function feat = gen_intron_features(gene, t, num_bins, read_len, reverse_ret)
% feat = gen_intron_features(gene, t, num_bins, read_len)
%
% -- input --
% gene: struct defining a gene with start, stops, exons etc.
% t: index of transcript
% num_bins: number of bins for PLiFs
% read_len: length of reads
%
% -- output --
% feat: PxDxD matrix of features at each position in the training set

myINF = 10^20; 

if nargin<5
  reverse_ret = 0;
end

offset = gene.start-1;
exons = gene.exons{t};

lmt = [ceil(linspace(0, sqrt(read_len), num_bins).^2) inf];

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

% indices relative to transcript and placed according to their exonic position 
gidx = zeros(length(eidx), 1);
gidx(tidx) = 1:length(tidx);

feat = zeros(length(eidx), num_bins, num_bins);
if size(exons,1)==1, % transcript with one exon
  feat(tidx, length(lmt)-1, length(lmt)-1) = feat(tidx, length(lmt)-1, length(lmt)-1) + 1;
else % transcript with multiple exons
  intron_bnd = exons(1:end-1,2)'-offset; % relative gene coordinates
  [tmp idx1 idx2] = intersect(intron_bnd, eidx);
  assert(isequal(intron_bnd, eidx(idx2)));
  intron_bnd = gidx(idx2) + 0.5; % relative to transcript and to their exonic position
  assert(all(intron_bnd>0));
  % distances to the closest intron start (in right or 3' direction)
  dist_starts = repmat(intron_bnd,1,length(tidx))-repmat(1:length(tidx),length(intron_bnd),1);
  dist_starts(dist_starts<0) = myINF;
  dist_starts = min(dist_starts, [], 1);
  % distances to the closest intron end (in left or 5' direction)
  dist_ends = repmat(1:length(tidx),length(intron_bnd),1)-repmat(intron_bnd,1,length(tidx));
  dist_ends(dist_ends<0) = myINF;
  dist_ends = min(dist_ends, [], 1);
  % reverse for minus strand
  if reverse_ret & gene.strand=='-'
    tmp = dist_starts; dist_starts = dist_ends; dist_ends = tmp;
    clear tmp;
  end
  for l = 1:length(lmt)-1,
    for r = 1:length(lmt)-1,
      fidx = find(lmt(l)<=dist_ends & lmt(l+1)>dist_ends & lmt(r)<=dist_starts & lmt(r+1)>dist_starts);
      feat(tidx(fidx), l, r) = feat(tidx(fidx), l, r) + 1;
    end
  end
end

% reverse for minus strand
if reverse_ret & gene.strand=='-'
  rev_idx = size(feat,1):-1:1;
  feat = feat(rev_idx,:,:);
end

assert(sum(sum(sum(feat,3),2),1)==tlen);
assert(all(sum(sum(feat,3),2)==1|sum(sum(feat,3),2)==0));
assert(size(feat,1)==length(eidx));