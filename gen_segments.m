function segments = gen_segments(gene)
% segments = gen_segments(gene)
%
% -- input --
% gene: struct defining a gene with start, stop, exons etc.
%
% -- output --
% segments: Sx2 vector of starts and stops in exonic coordinates
%           defining common segments S across exons of all transcripts

offset = gene.start-1;
eidx = gene.eidx;
starts = []; stops = [];
blocks = zeros(1,gene.stop-gene.start+1);
for t = 1:length(gene.transcripts),
  exons = gene.exons{t}-offset;    
  starts = [starts, exons(:,1)'];
  stops = [stops, exons(:,2)'];
  for e = 1:size(gene.exons{t},1),
    blocks(exons(e,1):exons(e,2)) = blocks(exons(e,1):exons(e,2)) + 1;
  end
end
starts = unique(starts); stops = unique(stops);
cand = stops(1:end-1)+1;
starts = unique([starts, cand(blocks(cand)>0)]);
cand = starts(2:end)-1;
stops = unique([stops, cand(blocks(cand)>0)]);

assert(all(starts<=stops));
assert(sum(stops-starts+1)==sum(blocks>0));

starts = starts + offset;
stops = stops + offset;

% in exonic coordinates
[tmp idx1 idx2] = intersect(starts, eidx);
assert(isequal(starts, eidx(idx2)));
starts = idx2;
[tmp idx1 idx2] = intersect(stops, eidx);
assert(isequal(stops, eidx(idx2)));
stops = idx2;

segments = [starts', stops'];