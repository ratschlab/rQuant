function [profile_weights, obj, seq_weights] = opt_density(CFG, genes)
% [profile_weights, obj, seq_weights] = opt_profiles_smo(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc.
%
% -- output --
% profile_weights: weights of profile functions 
% obj: objective evaluated with optimal parameters
% seq_weights: weights for sequence normalisation


T = length([genes.transcripts]);       % number of transcripts
P = sum([genes.exonic_len]);           % number of positions
P_all = P;
F = CFG.num_plifs;                     % number of supporting points
N = size(CFG.transcript_len_ranges,1); % number of length bins
if CFG.norm_seqbias
  S = 0;
  for o = 1:CFG.RR.order,
    S = S + (CFG.RR.half_win_size*2-o+1) * 4^o;
  end
end

I = 0; % upper bound for number of introns
%if norm_seqbias, Ps = 0; end % number of positions with sequence features
for g = 1:length(genes),
  introns = zeros(0,2);
  for t = 1:length(genes(g).transcripts),
    introns = [introns; genes(g).exons{t}(1:end-1,2)+1, genes(g).exons{t}(2:end,1)-1];
    if CFG.norm_seqbias,
      %Ps = Ps + sum(genes(g).exons{t}(:,2)-genes(g).exons{t}(:,1)+1) - 2*CFG.RR.half_win_size;
    end
  end
  introns = unique(introns, 'rows');
  I = I + size(introns,1);
  clear introns;
end

%%%%% pre-processing %%%%%
coverage = zeros(P, 1);
exon_feat = sparse(P, T); % stores exon features from all transcripts in profile gene set
exon_feat_val = zeros(P, 1);
exon_feat_row = [1:P]';
exon_feat_col = zeros(P, 1);
intron_count = zeros(I, 1);
intron_mask = zeros(I, T);
mask = true(P, 1); 
ci = 0; ct = 0; cn = 0;
if CFG.norm_seqbias
  cs = 0;
  seq_feat = sparse(P, S*T); % encoded sequence features
  %seq_target = sparse(1, P); % target values (read starts normalised by background)
  %seq_target_bg = sparse(1, P); % background (average number of read starts)
end
weights = zeros(1,T);
if CFG.VERBOSE>1, fprintf(1, 'Loading reads...\n'); tic; end
tmp_VERBOSE = CFG.VERBOSE;
CFG.VERBOSE = 0;
for g = 1:length(genes),
  fprintf('%i\r', g);
  try
    [tmp_coverage excluded_reads reads_ok tmp_introns] = get_coverage_per_read(CFG, genes(g), 1);
  catch
    reads_ok = 0;
  end
  assert(reads_ok==1);
  coverage(ci+[1:genes(g).exonic_len],1) = tmp_coverage;
  clear tmp_coverage;
  Tg = length(genes(g).transcripts);
  % initialisation of transcript weights to proportionate mean coverage
  weights(ct+[1:Tg]) = full(mean(coverage(ci+[1:genes(g).exonic_len]))/Tg*ones(1,Tg));
  for t = 1:length(genes(g).transcripts), % only works for tscp_len=1
    % exon features
    [exon_feat(ci+[1:genes(g).exonic_len], ct+1), tmp_exon_feat_val, feat_del_idx] = gen_exon_features(genes(g), t, F, CFG.max_side_len, 1);
    exon_feat_val(ci+[1:genes(g).exonic_len]) = tmp_exon_feat_val + F*ct; % transformation to linear indices
    exon_feat_col(ci+[1:genes(g).exonic_len]) = ct+1;
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
    clear tmp_exon_feat_val;
  end
  % sequence features
  if CFG.norm_seqbias
    genes(g).num_read_starts = tmp_read_starts;
    for t = 1:length(genes(g).transcripts),
      tmp_X = gen_sequence_features(CFG, genes(g), t);
      %[tmp_X tmp_Y tmp_Y_bg] = gen_sequence_features(CFG, genes(g), t);
      assert(size(tmp_X,1)==S);
      assert(size(tmp_X,2)<=genes(g).exonic_len);
      seq_feat(ci+[1:size(tmp_X,2)], cs+[1:S]) = tmp_X';
      cs = cs + S;
    end
    %seq_target(1, ci+[1:size(tmp_X,2)]) = tmp_Y;
    %seq_target_bg(1, ci+[1:size(tmp_X,2)]) = tmp_Y_bg;
    clear tmp_read_starts tmp_X tmp_Y tmp_Y_bg;
  end
  % introns
  [tmp_intron_mask tmp_intron_count] = get_intron_data(genes(g), CFG, tmp_introns, g);
  intron_mask(cn+[1:length(tmp_intron_count)], ct+[1:Tg]) = tmp_intron_mask; 
  intron_count(cn+[1:length(tmp_intron_count)],1) = tmp_intron_count;
  % repeat mask
  fname = sprintf('%s%s_repeat', CFG.repeats_fn, genes(g).chr);
  if exist(sprintf('%s.pos', fname), 'file')
    [map.pos map.repeats] = interval_query(fname, {'repeats'}, [genes(g).start;genes(g).stop]);
    if ~isempty(map.pos)
      [tmp idx1 idx2] = intersect(map.pos, genes(g).eidx);
      assert(length(idx2)<=length(map.pos));
      mask(ci+idx2) = false;
    end
  end
  mask(ci+feat_del_idx) = false;
  ci = ci + genes(g).exonic_len;
  ct = ct + length(genes(g).transcripts);
  cn = cn + length(tmp_intron_count);
  clear tmp_introns tmp_intron_mask tmp_intron_count feat_del_idx;
end
assert(P==size(exon_feat,1));
if CFG.norm_seqbias, assert(P==size(seq_feat,1)); end
% cut intron data to actual number of introns
if I < cn
  I = cn;
  assert(sum(sum(intron_mask(I+1:end,:)))==0);
  assert(sum(intron_count(I+1:end,1))==0);
  intron_mask = intron_mask(1:I,:);
  intron_count = intron_count(1:I,1);
end
CFG.VERBOSE = tmp_VERBOSE;
if CFG.VERBOSE>1, fprintf(1, 'Took %.1fs.\n', toc); end
tscp_len_bin = [genes.transcript_len_bin];
% find thetas that do not need to be optimised (located in the body of the profile function)
lmt = get_limits(CFG.max_side_len, F/2);
pw_nnz = false(F, N);
%max_tl = max([genes.transcript_length]);
%if max_tl<CFG.transcript_len_ranges(N,1)
%  CFG.transcript_len_ranges(N,2) = CFG.transcript_len_ranges(N,1) + 1;
%else
%  CFG.transcript_len_ranges(N,2) = max_tl;
%end
for n = 1:N,
  fidx = find(CFG.transcript_len_ranges(n,2)/2>lmt, 1, 'last');
  pw_nnz([1:fidx, F-fidx+1:F],n) = true;
end

% subsample positions
if CFG.subsample,
  subsample_frac = min(CFG.subsample_frac, CFG.max_num_train_exm/sum(mask));
  tmp_P = round(P*subsample_frac);
  midx = find(mask);
  ridx = randperm(length(midx));
  mask(midx(ridx(tmp_P+1:end))) = false;
  clear midx ridx;
end

% exclude positions that are repetitive or subsampled or not representable by PLiFs
if any(~mask),
  subs_idx = find(mask);
  P = length(subs_idx);
  coverage = coverage(subs_idx, :);
  exon_feat = exon_feat(subs_idx, :);
  exon_feat_val = exon_feat_val(subs_idx, 1);
  exon_feat_row = [1:P]';
  exon_feat_col = exon_feat_col(subs_idx, 1);
  if CFG.norm_seqbias
    seq_feat = seq_feat(subs_idx, :);
    %seq_target = seq_target(:, subs_idx);
    %seq_target_bg = seq_target_bg(:, subs_idx);
  end
  if CFG.VERBOSE>0, fprintf('subsampled from %i to %i positions\n', P_all, P); end
  clear P_old;
end

if 0
% exclude thetas that have too little example data
lmt(end) = max([genes.transcript_length]);
tmp1 = zeros(F,N);
for n = 1:N,
  for f = 1:F,
    idx_w_n = find(tscp_len_bin==n);
    idx_w_th = f+(idx_w_n-1)*F;
    tmp1(f,n) = sum(sum(exon_feat(:,idx_w_th)));
  end
end
diff = lmt(2:end)-lmt(1:end-1);
tmp2 = repmat([diff diff(end:-1:1)]', 1, N);
tmp3 = tmp1./tmp2;
pw_nnz(tmp3<(median(tmp3(tmp3>0))-std(tmp3(tmp3>0))) & tmp1>0) = false;
clear tmp1 tmp2 tmp3 diff;
else
pw_nnz = false(F, N);
for n = 1:N,
  fidx = find(CFG.transcript_len_ranges(n,2)/2>lmt, 1, 'last');
  pw_nnz([1:fidx+1, F-fidx:F],n) = true;
end
end
pw_nnz

out_dir = sprintf('%s/tmp/', CFG.out_dir);
if ~exist(out_dir ,'dir'),
  [s m mid] = mkdir(out_dir);
  assert(s);
end

%%%%% optimisation %%%%%
eps = 1e-2;
C_w = [genes.transcript_length]';
%%% initialisation of variables
weights_old = zeros(1,T);
if CFG.norm_seqbias
  seq_weights = ones(S,1);
  seq_weights_old = zeros(S, 1);
  num_opt_steps = 3;
else
  seq_weights = nan;
  num_opt_steps = 2;
end
profile_weights = zeros(F, N);
profile_weights(pw_nnz) = 1;
profile_weights_old = zeros(F, N);
pw_nnz = reshape(pw_nnz, 1, F*N);
fval = 1e100*ones(1, num_opt_steps); % objective values
fval_old = 0;
iter = 1;
tp_idx = zeros(1,F*T);
for n = 1:F,
  tp_idx([1:F:length(tp_idx)]+n-1) = n+[0:(F*T+F):F*T*T];
end
if CFG.norm_seqbias
  ts_idx = zeros(1,S*T);
  for n = 1:S,
    ts_idx([1:S:length(ts_idx)]+n-1) = n+[0:(S*T+S):S*T*T];
  end
end
lambda = 1;

if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\t\tDelta norm\n'); end
while 1
  tic
  weights_old = weights;
  fval_old = fval;
  profile_weights_old = profile_weights;
  if CFG.norm_seqbias
    seq_weights_old = seq_weights;
  end
  lambda_old = lambda;
  % pre-computations
  if CFG.norm_seqbias
    tmp_seq_weights = sparse(S*T,T);
    tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
    seq_coeff = seq_feat*tmp_seq_weights;
  else
    seq_coeff = [];
  end
  cnt = 1;
  
  %%%%% A. optimise transcript weights
  exon_mask = gen_exon_mask(profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_row, exon_feat_col);
  if CFG.norm_seqbias
    exon_mask = exon_mask.*seq_coeff;
  end
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  R_const = 1/lambda^2*(CFG.C_N*sum(sum((profile_weights(:,1:end-1)-profile_weights(:,2:end)).^2)) + CFG.C_F*sum(sum((profile_weights(1:end-1,:)-profile_weights(2:end,:)).^2)));
  [weights, fval(cnt)] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, lambda*intron_mask, lambda*C_w, R_const, 1, weights, 'L1');
  assert(fval_old(end)-fval(cnt)>-1e-3);
  cnt = cnt + 1;
  CFG.VERBOSE = tmp_VERBOSE;
    
  %%%%% B. optimise profile weights
  profile_weights = reshape(profile_weights, 1, F*N);
  R_const = abs(weights)*(lambda*C_w);
  if size(intron_mask,1)>0, R_const = R_const + CFG.C_I*sum((lambda*intron_mask*weights'-intron_count).^2); end
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  CFG.C_F = 1/lambda^2*CFG.C_F;
  CFG.C_N = 1/lambda^2*CFG.C_N;
  [profile_weights, fval(cnt)] = opt_profiles_descent(CFG, profile_weights, pw_nnz, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_row, exon_feat_col, weights, coverage, seq_coeff, R_const);
  assert(fval(cnt-1)-fval(cnt)>-1e-3);
  cnt = cnt + 1;
  CFG.VERBOSE = tmp_VERBOSE;
  
  %%%%% C. optimise sequence weights
  if CFG.norm_seqbias
    exon_mask = gen_exon_mask(profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_row, exon_feat_col);
    R_const = abs(weights)*(lambda*C_w) + 1/lambda^2*(CFG.C_N*sum(sum((profile_weights(:,1:end-1)-profile_weights(:,2:end)).^2)) + CFG.C_F*sum(sum((profile_weights(1:end-1,:)-profile_weights(2:end,:)).^2)));
    if size(intron_mask,1)>0, R_const = R_const + CFG.C_I*sum((lambda*intron_mask*weights'-intron_count).^2); end
    tmp_VERBOSE = CFG.VERBOSE;
    CFG.VERBOSE = 0;
    [seq_weights, fval(cnt)] = opt_seq_descent(CFG, seq_weights, seq_feat, exon_mask, tmp_seq_weights, weights, coverage, R_const);
    assert(fval(cnt-1)-fval(cnt)>-1e-3);
    cnt = cnt + 1;
    CFG.VERBOSE = tmp_VERBOSE;
    %plot_sequence_weights(seq_weights, CFG.RR.order, 2*CFG.RR.half_win_size);
  end
  
  if 0
  %%%% D. optimise regularisation parameters
  R1 = abs(weights)*C_w;
  R2 = CFG.C_N*sum(sum((profile_weights(:,1:end-1)-profile_weights(:,2:end)).^2)) + CFG.C_F*sum(sum((profile_weights(1:end-1,:)-profile_weights(2:end,:)).^2));
  if size(intron_mask,1)>0,
    R1 = R1 - 2*CFG.C_I*sum(intron_mask*weights'.*intron_count);
    R0 = CFG.C_I*sum((intron_mask*weights').^2);
    lambda = solve_lambda(R1, R2, R0);
  else
    lambda = solve_lambda(R1, R2);
  end
  if CFG.norm_seqbias
    tmp_seq_weights = sparse(S*T,T);
    tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
    seq_coeff = seq_feat*tmp_seq_weights;
    fval(cnt) = sum(((exon_feat*tmp_profiles).*seq_coeff*weights'-coverage).^2);
  else
    fval(cnt) = sum((exon_feat*tmp_profiles*weights'-coverage).^2);
  end
  fval(cnt) = fval(cnt) + abs(weights)*(lambda*C_w) + 1/lambda^2*(CFG.C_N*sum(sum((profile_weights(:,1:end-1)-profile_weights(:,2:end)).^2)) + CFG.C_F*sum(sum((profile_weights(1:end-1,:)-profile_weights(2:end,:)).^2)));
  if size(intron_mask,1)>0, fval(cnt) = fval(cnt) + CFG.C_I*sum((lambda*intron_mask*weights'-intron_count).^2); end
  end
  
  %%%%% convergence criteria
  if CFG.norm_seqbias
    norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N), seq_weights_old', lambda_old] - [weights, reshape(profile_weights,1,F*N), seq_weights' lambda]);
  else
    norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N), lambda_old] - [weights, reshape(profile_weights,1,F*N), lambda]);
  end
  if fval_old(end)>=fval(end), sg = '-'; else sg = '+'; end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.3d\t%.5f\t%.3d\t%s\t%.1f\n', iter, fval(end), fval(end)/fval_old(end), norm_weights, sg, toc); end
  
  if norm(fval_old-fval)<eps || norm_weights<eps || iter >= CFG.max_iter
    break;
  end
  %profile_weights
  if iter>15, keyboard; end
  save(sprintf('%siter%i.mat', out_dir, iter), 'profile_weights', 'weights', 'lambda', 'seq_weights');
  iter = iter + 1;
end
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

% normalise weights
if 0
for n = 1:N,
  profile_weights(:,n) = profile_weights(:,n)./norm(profile_weights(:,n));
end
seq_weights  = seq_weights./norm(seq_weights);
end

obj = fval(end);
