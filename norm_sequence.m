function coverage_norm = norm_sequence(CFG, gene, coverage)
% coverage_norm = norm_sequence(CFG, gene, coverage)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc. 
% coverage: matrix of exonic positions x reads
%
% -- output --
% coverage_norm: matrix of exonic positions x reads, normalised 

fidx = find((gene.eidx(2:end)-gene.eidx(1:end-1))>1);
exons = zeros(length(fidx)+1, 2);
exons(:,1) = [gene.eidx(1); gene.eidx(fidx+1)'];
exons(:,2) = [gene.eidx(fidx)'; gene.eidx(end)];
seq = load_genomic(gene.chr, gene.strand, exons(:,1), exons(:,2), CFG.genome_info, 0);
assert(length(seq)==length(gene.num_read_starts));
% input sequence data for regression
X = char(zeros(CFG.RR.half_win_size*2, length(seq)-2*CFG.RR.half_win_size));
for x = CFG.RR.half_win_size+1:length(seq)-CFG.RR.half_win_size,
  X(:,x-CFG.RR.half_win_size) = upper(char(seq(x-CFG.RR.half_win_size:x+CFG.RR.half_win_size-1)));
end
if ~isnan(CFG.RR.seq_norm_weights)
  Y_pred = predict_Ridge(CFG, X, CFG.RR.seq_norm_weights);
  num_read_starts_pred = [gene.num_read_starts(1:CFG.RR.half_win_size), exp(Y_pred) * sum(gene.num_read_starts) / length(gene.num_read_starts), gene.num_read_starts(end-CFG.RR.half_win_size:end)]; 
else
  num_read_starts_pred = gene.num_read_starts;
end
diff = coverage - [zeros(1,size(coverage,2)); coverage(1:end-1,:)];
[m midx] = max(diff, [], 1);
% weights each read according to the corrected read start frequency
coverage_norm = double(coverage) .* repmat(num_read_starts_pred(midx)./gene.num_read_starts(midx), size(coverage,1), 1);