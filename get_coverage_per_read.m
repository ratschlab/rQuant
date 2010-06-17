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

if nargin<4
  reverse_ret = 0;
end

if ~isfield(CFG, 'both_strands')
  CFG.both_strands = 0;
end

win = CFG.read_len;

eidx = [max(gene.eidx(1)-win,1):gene.eidx(1)-1, gene.eidx, gene.eidx(end)+1:min(gene.eidx(end)+win,CFG.chr_len(gene.chr_num))];
win_size = length(max(gene.eidx(1)-win,1):gene.eidx(1)-1);

for f = 1:length(CFG.read_maps_fn{gene.chr_num}),
  fname = CFG.read_maps_fn{gene.chr_num}{f};
  if ~fexist(fname),
    warning('BAM file %s does not exist', fname);
  end
  try
    if CFG.both_strands
      [coverage_idx_tmp{f}, intron_list_tmp] = get_reads(fname, gene.chr, eidx(1), eidx(end), '0');
    else
      [coverage_idx_tmp{f}, intron_list_tmp] = get_reads(fname, gene.chr, eidx(1), eidx(end), gene.strand); 
    end
  catch
    warning('get_reads failed');
    intron_list = intron_list';
    ok = 0;
    return;
  end
  if ~isempty(intron_list_tmp),
    intron_list = [intron_list intron_list_tmp{:}];
  end
end

% process coverage: convert to exonic position indices
coverage_idx = [coverage_idx_tmp{:}];
if ~isempty(coverage_idx)
  coverage = sparse(coverage_idx(1,:), coverage_idx(2,:), 1, max(coverage_idx(1,:)), eidx(end)-eidx(1)+1)';
else
  coverage = sparse([], [], 1, eidx(end)-eidx(1)+1, 0);
end
coverage = coverage(eidx(win_size+1:win_size+gene.exonic_len)-eidx(win_size+1)+1, :);
% no overlapping reads
assert(~any(any(full(coverage>1))));

% process intron list (1: intron start, 2: intron stop, 3: confirmation, 4: strand)
intron_list = [intron_list', zeros(size(intron_list,2), 1), (gene.strand=='-')*ones(size(intron_list,2), 1)];
intron_list_unique = unique(intron_list, 'rows');
for n = 1:size(intron_list_unique,1),
  intron_list_unique(n,3) = sum(intron_list_unique(n,1)==intron_list(:,1) & ...
                                intron_list_unique(n,2)==intron_list(:,2) & ...
                                intron_list_unique(n,4)==intron_list(:,4));
end
intron_list = intron_list_unique;
clear intron_list_unique;

% process read starts: count reads starting at each exonic position
for c = 1:size(coverage, 2),
  fidx = find(coverage(:,c)~=0, 1, 'first');
  if ~isempty(fidx),
    read_starts(fidx) = read_starts(fidx) + 1;
  end
end

% determine set of excluded reads
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

% reverse for minus strand
if reverse_ret && gene.strand=='-'
  rev_idx = size(coverage,1):-1:1;
  coverage = coverage(rev_idx,:);
end