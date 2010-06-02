function [mask excluded_reads ok intron_list read_starts_pos] = get_coverage_per_read(CFG, gene, reverse_ret)
% [mask excluded_reads ok intron_list read_starts] = get_coverage_per_read(CFG, gene, reverse_ret)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc.
% reverse_ret: reverse exon_count for minus strand
%
% -- output --
% mask: matrix of exonic positions x reads
% excluded_reads: reads exluded per transcript
% intron_list: nx4 list of introns (intron start, intron stop, confirmation, strand)
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
intron_list = zeros(2,0);
read_starts = zeros(1,0);
read_starts_pos = [];
for f = 1:length(CFG.read_maps_fn{gene.chr_num}),
  fname = CFG.read_maps_fn{gene.chr_num}{f};
  if ~fexist(fname),
    warning('BAM file %s does not exist', fname);
  end
  try
    if CFG.both_strands
      [mask_idx, read_intron_list] = get_reads(fname, gene.chr, eidx(1),eidx(end),'0');
    else
      [mask_idx, read_intron_list] = get_reads(fname, gene.chr, eidx(1),eidx(end),gene.strand); 
    end
    if ~isempty(mask_idx)
      mask_tmp{f} = sparse(mask_idx(1,:), mask_idx(2,:), 1, max(mask_idx(1,:)),eidx(end)-eidx(1)+1);
      mask_tmp{f} = full(mask_tmp{f}(:,eidx-eidx(1)+1)); % remove this part to continue working with a sparse matrix
    else
      mask_tmp{f} = zeros(0,length(eidx));
    end
  catch
    warning('get_reads failed');
    intron_list = intron_list';
    ok = 0;
    return;
  end
  
  if ~isempty(mask_tmp{f}),
    intron_list = [intron_list read_intron_list{:}];
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
  read_starts_ = zeros(1, size(mask_tmp{f},1));
  for i=1:size(mask_tmp{f}, 1),
    idx = find(mask_tmp{f}(i,:)~=0, 1, 'first');
    if ~isempty(idx),
      read_starts_(i) = idx;
    else
      read_starts_(i) = nan;
    end
  end
  read_starts = [read_starts read_starts_(~isnan(read_starts_))];
end

% process intron list
intron_list = intron_list';
intron_list_unique = unique(intron_list, 'rows');
intron_list(:,3)=0;
intron_list(:,4)=0;
intron_list_unique(:,3)=0;
intron_list_unique(:,4)=0;
for i=1:size(intron_list_unique,1),
  cnt = sum(intron_list_unique(i,1)==intron_list(:,1) & ...
            intron_list_unique(i,2)==intron_list(:,2) & ...
            intron_list_unique(i,4)==intron_list(:,4));
  intron_list_unique(i,3)=cnt;
end
intron_list=intron_list_unique;
clear intron_list_unique;

% process read_start list
read_starts_pos = zeros(1, gene.exonic_len);
for i=1:length(read_starts),
  read_starts_pos(read_starts(i)) = read_starts_pos(read_starts(i)) + 1;
end

% collect masks
%mask = zeros(gene.exonic_len, num_reads);
mask = false(gene.exonic_len, num_reads);
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
eidx = 1:length(gene.eidx);
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