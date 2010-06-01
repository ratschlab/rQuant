function profile_norm = norm_sequence(CFG, gene, transcript_idx, profile)
% coverage_norm = norm_sequence(CFG, gene, coverage)
%
% -- input --
% CFG: configuration struct
% gene: struct defining a gene with start, stops, exons etc. 
%
% -- output --
% coverage_norm: matrix of exonic positions x reads, normalised 

exons = gene.exons{transcript_idx};
exons(1,1) = exons(1,1)-CFG.RR.half_win_size; 
exons(end,2) = exons(end,2) + CFG.RR.half_win_size + 1;

if ~isfield(gene, 'strands') || length(gene.strands)<transcript_idx,
  gene.strands(transcript_idx) = gene.strands;
end
seq = load_genomic(gene.chr, gene.strands(transcript_idx), exons(:,1), exons(:,2), CFG.genome_info, 0);
assert(length(seq)==length(profile)+2*CFG.RR.half_win_size+1);

% input sequence data for regression
X = char(zeros(CFG.RR.half_win_size*2, length(seq)-2*CFG.RR.half_win_size));
for x = CFG.RR.half_win_size+1:length(seq)-CFG.RR.half_win_size,
  X(:, x-CFG.RR.half_win_size) = upper(char(seq(x-CFG.RR.half_win_size:x+CFG.RR.half_win_size-1)));
end
if ~any(isnan(CFG.RR.seq_norm_weights)),
  Y_pred = predict_Ridge(CFG, X, CFG.RR.seq_norm_weights);
  num_read_starts_pred = exp(Y_pred);
else
  num_read_starts_pred = ones(1, length(seq)-2*CFG.RR.half_win_size-1);
end

coverage = zeros(1, length(num_read_starts_pred));
for i=1:length(num_read_starts_pred),
  idx = i:min(i+CFG.read_len-1, length(num_read_starts_pred));
  coverage(idx) = coverage(idx) + num_read_starts_pred(i);
end
coverage = coverage / mean(coverage);

profile_norm = profile .* coverage;

