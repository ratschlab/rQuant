function w = train_norm_sequence(CFG, genes)
% seq_norm_weights = train_norm_sequence(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc. 
%
% -- output --
% w: weight vector of trained Ridge regression


w = nan; 

if CFG.VERBOSE>0, fprintf(1, 'Generating training set for sequence normalisation...\n'); tic; end
TR.X = []; TR.Y = [];  
for g = 1:length(genes),  
  gene = genes(g);
  exons = gene.exons;
  eidx = [];
  for e = 1:size(exons,1),
    eidx = [eidx, exons(e,1):exons(e,2)];
  end
  eidx = unique(eidx);
  [tmp idx1 idx2] = intersect(eidx, gene.eidx);
  num_read_starts = gene.num_read_starts(idx2);
  seq = load_genomic(gene.chr, gene.strand, exons(:,1), exons(:,2), CFG.genome_info, 0);
  assert(length(seq)==length(num_read_starts));

  num_read_starts_norm = length(num_read_starts) * num_read_starts / sum(num_read_starts) + 1e-2;
  num_read_starts_norm(1:CFG.RR.half_win_size) = nan;
  num_read_starts_norm(end-CFG.RR.half_win_size:end) = nan;
  idx = find(~isnan(num_read_starts_norm));

  % target values for regression
  TR.Y = [TR.Y log(num_read_starts_norm(idx))'];
  % input sequence data for regression
  X = char(zeros(CFG.RR.half_win_size*2, length(idx)));
  for x = 1:length(idx),
    X(:,x) = char(seq(idx(x)-CFG.RR.half_win_size:idx(x)+CFG.RR.half_win_size-1));
  end
  TR.X = [TR.X upper(X)];
end
if CFG.VERBOSE>0, fprintf(1, 'Took %.1fs.\n', toc); end

ridx = randperm(size(TR.X,2));
num_train = floor(length(ridx)*CFG.RR.num_train_frac);

if num_train < 10^4
  if CFG.VERBOSE>0, fprintf(1, 'insufficient number of training examples (%i) for Ridge regression\n', num_train); end
  return;
end

X_train = TR.X(:,ridx(1:num_train));
Y_train = TR.Y(:,ridx(1:num_train));
X_eval = TR.X(:,ridx(num_train+1:end));
Y_eval = TR.Y(:,ridx(num_train+1:end));

% train weights of Ridge regression
if CFG.VERBOSE>0, fprintf(1, '%i training examples for Ridge regression\n', size(X_train,2)); end
[w TR] = train_Ridge(CFG, X_train, Y_train);
% test on omitted examples
if CFG.VERBOSE>0,
  fprintf(1, '%i test examples for Ridge regression\n', size(X_eval,2));
  Y_pred = predict_Ridge(CFG, X_eval, w);
  % absolute variablity on test set
  TE.Q1 = mean(abs(Y_pred - Y_eval)) / mean(abs(Y_eval - median(Y_eval)));
  % squared variablity on test set
  TE.Q2 = mean((Y_pred - Y_eval).^2) / mean((Y_eval - mean(Y_eval)).^2);
  fprintf(1, 'relative variability (Q1): training %.3f, test %.3f\n', TR.Q1, TE.Q1);
  fprintf(1, 'relative variability (Q2): training %.3f, test %.3f\n', TR.Q2, TE.Q2);
end