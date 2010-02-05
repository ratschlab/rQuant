function [mask excluded_reads ok] = get_coverage_per_read(CFG, gene, reverse_ret)
% [mask excluded_reads ok] = get_coverage_per_read(CFG, gene, reverse_ret)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc.
% reverse_ret: reverse exon_count for minus strand
%
% -- output --
% mask: matrix of exonic positions x reads
% excluded_reads: reads exluded per transcript
% ok: indicates success of file parsing

MIN_COV = 0.5;
ok = 1;
mask = [];
excluded_reads = [];

if nargin<4
  reverse_ret = 0;
end

if ~isfield(CFG, 'both_strands')
  CFG.both_strands = 0;
end

win = CFG.read_len;

eidx = [max(gene.eidx(1)-win,1):gene.eidx(1)-1, gene.eidx, gene.eidx(end)+1:min(gene.eidx(end)+win,CFG.chr_len(gene.chr_num))];
win_size = length(max(gene.eidx(1)-win,1):gene.eidx(1)-1);

num_reads = 0;
file_exist = zeros(1,length(CFG.read_maps_fn{gene.chr_num}));
for f = 1:length(CFG.read_maps_fn{gene.chr_num}),
  fname = CFG.read_maps_fn{gene.chr_num}{f};
  try
    if CFG.both_strands
      mask_tmp{f} = get_reads(fname, CFG.samtools_dir, gene.chr, '0', eidx, CFG.paired);
    else
      mask_tmp{f} = get_reads(fname, CFG.samtools_dir, gene.chr, gene.strand, eidx, CFG.paired); 
    end
    mask_tmp{f}(mask_tmp{f}>1) = 1;
    % delete rows with zeros and those with overlapping reads
    mask_tmp{f}(sum(mask_tmp{f},2)==0 | any(mask_tmp{f}>1,2),:) = [];
    mask_tmp{f} = mask_tmp{f}(:,win_size+1:win_size+gene.exonic_len);
    assert(all(all(mask_tmp{f}==0 | mask_tmp{f}==1)));
    num_reads = num_reads + size(mask_tmp{f},1);
    if ~isempty(mask_tmp{f})
      file_exist(f) = 1;
    end
  catch
    ok = 0;
    return;
  end
end

% collect masks
%mask = zeros(gene.exonic_len, num_reads);
mask = logical(zeros(gene.exonic_len, num_reads));
cnt = 0;
for f = 1:length(CFG.read_maps_fn{gene.chr_num}),
  if file_exist(f)
    mask(:,cnt+1:cnt+size(mask_tmp{f},1)) = mask_tmp{f}';
    cnt = cnt + size(mask_tmp{f},1);
  end
end
assert(num_reads==cnt);
clear mask_tmp;

offset = gene.start-1;
eidx = [1:length(gene.eidx)];
sum_all = sum(mask,1);
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
  excluded_reads{t} = find(sum(mask(tidx,:),1)./sum_all<MIN_COV);
end

% reverse for minus strand
if reverse_ret && gene.strand=='-'
  rev_idx = size(mask,1):-1:1;
  mask = mask(rev_idx,:);
end