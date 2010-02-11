function w = train_norm_sequence(CFG, genes)
% seq_norm_weights = train_norm_sequence(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc. 
%
% -- output --
% w: weight vector of trained Ridge regression


num_nt = 4; num_feat = 0;
for o = 1:CFG.RR.order,
  num_feat = num_feat + (2*CFG.RR.half_win_size-o+1) * num_nt^o;
end

w = ones(num_feat,1); 

if CFG.VERBOSE>0, fprintf(1, 'Generating training set for sequence normalisation...\n'); tic; end
TR.X = []; TR.Y = [];  
for g = 1:length(genes),  
  gene = genes(g);
  exons = gene.exons;
  eidx = [];
  for e = 1:size(exons,1),
    eidx = [eidx, exons(e,1):exons(e,2)];
    eidx = unique(eidx);
  end
  [tmp idx1 idx2] = intersect(eidx, gene.eidx);
  num_read_starts = gene.num_read_starts(idx2);
  seq = load_genomic(gene.chr, gene.strand, exons(:,1), exons(:,2), CFG.genome_info, 0);
  assert(length(seq)==length(num_read_starts));
  num_read_starts = length(num_read_starts)*num_read_starts/sum(num_read_starts)+1e-2;
  num_read_starts(1:CFG.RR.half_win_size) = nan;
  num_read_starts(end-CFG.RR.half_win_size:end) = nan;
  idx = find(~isnan(num_read_starts));
  TR.Y = [TR.Y log(num_read_starts(idx))];
  X = char(zeros(CFG.RR.half_win_size*2, length(idx)));
  for j = 1:length(idx),
    X(:,j) = char(seq(idx(j)-CFG.RR.half_win_size:idx(j)+CFG.RR.half_win_size-1));
  end
  TR.X = [TR.X upper(X)];
end
if CFG.VERBOSE>0, fprintf(1, 'Took %.1fs.\n', toc); end

ridx = randperm(size(TR.X,2));
num_train = floor(length(ridx)*CFG.RR.num_train_frac);

if num_train < 10^6
  if CFG.VERBOSE>0, fprintf(1, 'insufficient number of training examples (%i) for Ridge regression\n', num_train); end
  return;
end

X_train = TR.X(:,ridx(1:num_train));
Y_train = TR.Y(:,ridx(1:num_train));
X_eval = TR.X(:,ridx(num_train+1:end));
Y_eval = TR.Y(:,ridx(num_train+1:end));

[w TR] = train_Ridge(CFG, X_train, Y_train);
Y_pred = predict_Ridge(CFG, X_eval, w);
TE.Q1 = mean(abs(Y_pred - Y_eval)) / mean(abs(Y_eval - median(Y_eval)));
TE.Q2 = mean((Y_pred - Y_eval).^2) / mean((Y_eval - mean(Y_eval)).^2);

if CFG.VERBOSE>0,
  fprintf(1, 'absolute variability: training %.3f, test %.3f\n', TR.Q1, TE.Q1);
  fprintf(1, 'squared variability: training %.3f, test %.3f\n', TR.Q2, TE.Q2);
end