function [mask excluded_reads ok] = get_coverage_per_read(CFG, gene, flag, reverse_ret)
% [mask excluded_reads ok] = get_coverage_per_read(CFG, gene, flag, reverse_ret)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc.
% flag: enables single read parsing
% reverse_ret: reverse exon_count for minus strand
%
% -- output --
% mask: logical matrix of exonic positions x reads
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

if CFG.both_strands
  file_exist = zeros(1,3);
else
  file_exist = zeros(1,2);
end

% mapped reads
fname = sprintf('%s%s+_mapped.bam', CFG.read_maps_fn, gene.chr);
if exist(fname, 'file')
  try
    mask1 = get_reads(fname, CFG.samtools_dir, gene.chr, '+', eidx, flag);
    mask1(mask1>1) = 1;
    % delete rows with zeros and those with overlapping reads
    mask1(sum(mask1,2)==0 | any(mask1>1,2),:) = [];
    mask1 = mask1(:,win_size+1:win_size+gene.exonic_len);
    assert(all(all(mask1==0 | mask1==1)));
    num_reads = num_reads + size(mask1,1);
    if ~isempty(mask1)
      file_exist(1) = 1;
    end
  catch
    ok = 0;
    return;
  end
end
% spliced reads (strand-specific)
fname = sprintf('%s%s+_%s_spliced.bam', CFG.read_maps_fn, gene.chr, CFG.read_maps_select);
if exist(fname, 'file')
  try
    mask2 = get_reads(fname, CFG.samtools_dir, gene.chr, gene.strand, eidx, flag);
    mask2(mask2>1) = 1;
    % delete rows with zeros and those with overlapping reads
    mask2(sum(mask2,2)==0 | any(mask2>1,2),:) = [];
    mask2 = mask2(:,win_size+1:win_size+gene.exonic_len);
    assert(all(all(mask2==0 | mask2==1)));
    num_reads = num_reads + size(mask2,1);
    if ~isempty(mask2)
      file_exist(2) = 1;
    end
  catch
    ok = 0;
    return 
  end
end
% spliced reads (other strand)  
if CFG.both_strands
  if gene.strand=='+'
    other_strand = '-';
  else
    other_strand = '+';
  end
  fname = sprintf('%s%s+_%s_spliced.bam', CFG.read_maps_fn, gene.chr, CFG.read_maps_select);
  if exist(fname, 'file')
    try
      mask3 = get_reads(fname, CFG.samtools_dir, gene.chr, other_strand, eidx, flag);
      mask3(mask3>1) = 1;
      % delete rows with zeros and those with overlapping reads
      mask3(sum(mask3,2)==0 | any(mask3>1,2),:) = [];
      mask3 = mask3(:,win_size+1:win_size+gene.exonic_len);
      assert(all(all(mask3==0 | mask3==1)));
      num_reads = num_reads + size(mask3,1);
      if ~isempty(mask3)
        file_exist(3) = 1;
      end
    catch
      ok = 0;
      return;
    end
  end
end
% collect masks
mask = logical(zeros(gene.exonic_len, num_reads));
cnt = 0;
if file_exist(1)
  mask(:,cnt+1:cnt+size(mask1,1)) = mask1';
  cnt = cnt + size(mask1,1);
end
if file_exist(2)
  mask(:,cnt+1:cnt+size(mask2,1)) = mask2';
  cnt = cnt + size(mask2,1);
end
if CFG.both_strands && file_exist(3)
  mask(:,cnt+1:cnt+size(mask3,1)) = mask3';
  cnt = cnt + size(mask3,1);
end
assert(num_reads==cnt);
clear mask1 mask2 mask3;

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