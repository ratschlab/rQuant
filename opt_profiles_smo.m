function [profile_weights, obj, seq_norm_weights] = opt_profiles_smo(CFG, genes)
% [profile_weights, obj, seq_norm_weights] = opt_profiles_smo(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc.
%
% -- output --
% profile_weights: weights of profile functions 
% obj: objective evaluated with optimal parameters
% seq_norm_weights: weight vector of trained Ridge regression for sequence normalisation


T = length([genes.transcripts]);       % number of transcripts
P = sum([genes.exonic_len]);           % number of positions
P_all = P;
F = CFG.num_plifs;                     % number of supporting points
N = size(CFG.transcript_len_ranges,1); % number of length bins
if CFG.norm_seqbias
  S = 0;
  for o = 1:CFG.RR.order,
    S = S + (CFG.RR.half_win_size*2-o+1) * 4^o;;
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
end

%%%%% pre-processing
coverage = sparse(P, 1);
exon_feat = sparse(P, F*T); % stores exon features from all transcripts in profile gene set
intron_count = zeros(I, 1);
intron_mask = zeros(I, T);
mask = true(P, 1); 
ci = 0; cj = 0; ct = 0; cn = 0;
if CFG.norm_seqbias
  seq_feat = sparse(S, P);
  seq_target = sparse(1, P);
end
weights = zeros(1,T);
if CFG.VERBOSE>1, fprintf(1, 'Loading reads...\n'); tic; end
tmp_VERBOSE = CFG.VERBOSE;
CFG.VERBOSE = 0;
for g = 1:length(genes),
  try
    if CFG.norm_seqbias
      [tmp_coverage excluded_reads reads_ok tmp_introns tmp_read_starts] = get_coverage_per_read(CFG, genes(g), 1);
    else
      [tmp_coverage excluded_reads reads_ok tmp_introns] = get_coverage_per_read(CFG, genes(g), 1);
    end
  catch
    reads_ok = 0;
  end
  assert(reads_ok==1);
  coverage(ci+[1:genes(g).exonic_len],1) = tmp_coverage;
  clear tmp_coverage;
  Tg = length(genes(g).transcripts);
  % initialisation of transcript weights to proportionate mean coverage
  weights(ct+[1:Tg]) = full(mean(coverage(ci+[1:genes(g).exonic_len]))/Tg*ones(1,Tg));
  for t = 1:length(genes(g).transcripts),
    % exon features
    exon_feat(ci+[1:genes(g).exonic_len], cj+[1:F]) = gen_exon_features(genes(g), t, F, CFG.max_side_len, 1);
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
    cj = cj + F;
  end
  % sequence features
  if CFG.norm_seqbias
    genes(g).num_read_starts = tmp_read_starts;
    [tmp_X tmp_Y] = gen_sequence_features(CFG, genes(g), t);
    assert(size(tmp_X,1)==S);
    assert(size(tmp_X,2)==genes(g).exonic_len);
    seq_feat(:, ci+[1:size(tmp_X,2)]) = tmp_X;
    seq_target(1, ci+[1:size(tmp_X,2)]) = tmp_Y;
    clear tmp_read_starts tmp_X tmp_Y;
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
  ci = ci + genes(g).exonic_len;
  ct = ct + length(genes(g).transcripts);
  cn = cn + length(tmp_intron_count);
  clear tmp_introns tmp_intron_mask tmp_intron_count;
end
assert(P==size(exon_feat,1));
if CFG.norm_seqbias, assert(P==size(seq_feat,2)); end
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
lmt = linspace(0, sqrt(CFG.max_side_len), (F/2)+1).^2;
lmt(end) = inf;
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

% exclude positions that are repetitive or subsampled
if any(~mask),
  subs_idx = find(mask);
  P = length(subs_idx);
  coverage = coverage(subs_idx, :);
  exon_feat = exon_feat(subs_idx, :);
  if CFG.norm_seqbias
    seq_feat = seq_feat(:, subs_idx);
    seq_target = seq_target(:, subs_idx);
  end
  fprintf('subsampled from %i to %i positions\n', P_all, P);
  clear P_old;
end

% TODO
%   subsampling of positions (also repeat positions)
%   speed up
%%%%% optimisation
eps = 1e-3;
C_w = [genes.transcript_length]';
% initialisation of variables
weights_old = zeros(1,T);
profile_weights = zeros(F, N);
profile_weights(pw_nnz) = (F*N)/sum(sum(pw_nnz));
pw_nnz = reshape(pw_nnz, 1, F*N);
profile_weights_old = zeros(F, N);
if CFG.norm_seqbias
  seq_norm_weights = ones(S, 1);
  seq_norm_weights_old = ones(S, 1);
else
  seq_norm_weights = nan;
end
fval = 1e100; %1e100*ones(1,T+F*N);
fval_old = 0;
iter = 1;
tp_idx = zeros(1,F*T);
for n = 1:F,
  tp_idx([1:F:length(tp_idx)]+n-1) = n+[0:(F*T+F):F*T*T];
end
if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tDelta norm\n'); end
while 1
  tic
  weights_old = weights;
  fval_old = fval;
  profile_weights_old = profile_weights;
  if CFG.norm_seqbias
    seq_norm_weights_old = seq_norm_weights;
  end
  % pre-computations
  tmp_profiles = sparse(zeros(F*T,T));
  tmp_profiles(tp_idx) = profile_weights(:,tscp_len_bin);
  
  % A. optimise transcript weights
  exon_mask = exon_feat*tmp_profiles;
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  [weights, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, 1, weights, 'L1');
  %corr(weights', [genes.expr_orig]')
  CFG.VERBOSE = tmp_VERBOSE;
    
  % B. optimise profile weights
  profile_weights = reshape(profile_weights, 1, F*N);
  num_changed = 0; examine_all = true;
  changed = zeros(1, F*N);
  pidx = find(pw_nnz(1:N*F-1));
  ii = 0;
  for p = pidx,
    qidx = find(pw_nnz);
    qidx(qidx<=p) = [];
    ridx = randperm(length(qidx));
    qidx = qidx(ridx);
    for q = qidx,
      ii = ii + 1;
      %fprintf('%3.2f\r', 100*ii/sum([1:sum(pw_nnz)-1]));
      theta1 = profile_weights(p);
      theta2 = profile_weights(q);
      d = theta1 + theta2;
      f1 = mod(p,F); if f1==0, f1 = F; end
      f2 = mod(q,F); if f2==0, f2 = F; end
      n1 = ceil(p/F);
      n2 = ceil((q)/F);
      if n1~=n2, continue; end % theta pair is from same transcript length bin
      idx_w_n1 = find(tscp_len_bin==n1);  % all transcripts in length bin n1
      idx_w_n2 = find(tscp_len_bin==n2);  % all transcripts in length bin n2
      idx_w_th1 = f1+(idx_w_n1-1)*F;      % entries corresponding to theta1
      idx_w_th2 = f2+(idx_w_n2-1)*F;      % entries corresponding to theta2
      idx_wo_n1 = find(tscp_len_bin~=n1); % all transcripts not in n1
      idx_wo_n2 = find(tscp_len_bin~=n2); % all transcripts not in n2
      idx_wo_th1 = f1+(idx_wo_n1-1)*F;    % entries corresponding to theta_f1_notn1
      idx_wo_th2 = f2+(idx_wo_n2-1)*F;    % entries corresponding to theta_f2_notn2
      idx_wo_th = setdiff([1:F*T], [f1:F:F*T, f2:F:F*T]);
      % residue for theta2
      Rth2 = exon_feat(:,idx_w_th2)*weights(idx_w_n2)' - exon_feat(:,idx_w_th1)*weights(idx_w_n1)';
      % residue for theta1/theta2 independent variables
      R1 = exon_feat(:,idx_wo_th1)*tmp_profiles(idx_wo_th1,idx_wo_n1)*weights(idx_wo_n1)' + ...
           exon_feat(:,idx_wo_th2)*tmp_profiles(idx_wo_th2,idx_wo_n2)*weights(idx_wo_n2)' + ...
           exon_feat(:,idx_wo_th)*tmp_profiles(idx_wo_th,:)*weights' + ...
           exon_feat(:,idx_w_th1)*weights(idx_w_n1)'*d - ...
           coverage;
      R2 = 0;
      if n1<N, R2 = R2 + (profile_weights(f1+n1*F)-d)^2; end
      if n2<N, R2 = R2 + profile_weights(f2+n2*F)^2; end
      for f = 1:F,
        for n = 1:N-1,
          if n~=n1 & n~=n2, R2 = R2 + (profile_weights(f+(n-1)*F)-profile_weights(f+n*F))^2; end
        end
      end
      R3 = 0;
      if f1<F, R3 = R3 + (profile_weights(f1+1+(n1-1)*F)-d)^2; end
      if f2<F, R3 = R3 + profile_weights(f2+1+(n2-1)*F)^2; end
      for n = 1:N,
        for f = 1:F-1,
          if f~=f1 & f~=f2, R2 = R2 + (profile_weights(f+(n-1)*F)-profile_weights(f+1+(n-1)*F))^2; end
        end
      end  
      S1 = 2*sum(Rth2.^2) + 8; % quadratic term
      assert(S1>0); % condition for minimum (2nd derivative > 0)
      S2 = 4*d - 2*sum(Rth2'*R1); % linear term
      % coupling constraints
      if n1<N, S2 = S2 - 2*profile_weights(f1+n1*F); end % theta_f1,n1+1
      if n2<N, S2 = S2 + 2*profile_weights(f2+n2*F); end % theta_f2,n2+1
      if f1<F, S2 = S2 - 2*profile_weights(f1+1+(n1-1)*F); end % theta_f1+1,n1
      if f2<F, S2 = S2 + 2*profile_weights(f2+1+(n2-1)*F); end % theta_f2+1,n2
      S3 = sum(R1.^2) + R2 + R3; % constant term
      theta2_new = S2/S1;
      % clipping of theta2
      if theta2_new<0
        theta2_new = 0.0;
      end
      if theta2_new>d
        theta2_new = d;
      end
      theta1_new = d - theta2_new;
      % check if thetas have been changed
      if abs(theta1_new-theta1)>eps
        num_changed = num_changed + 1;
        changed(p) = changed(p) + 1;
      end 
      if abs(theta2_new-theta2)>eps
        num_changed = num_changed + 1;
        changed(q) = changed(q) + 1;
      end
      % update thetas and objective value
      profile_weights(p) = theta1_new;
      profile_weights(q) = theta2_new;
      fval(end+1) = quad_fun(theta2_new, S1/2, -S2, S3);
      tmp_pw = reshape(profile_weights, F, N);
      tmp_profiles(tp_idx) = tmp_pw(:,tscp_len_bin);
    end
  end
  profile_weights = reshape(profile_weights, F, N)
  %figure(iter); plot(profile_weights); ylim([0 10]);
  %plot(iter, fval(end), 'x');
  
  % C. determine sequence bias
  if CFG.norm_seqbias
    exon_mask = exon_feat*tmp_profiles;
    seq_target = log(exp(seq_target)./(sum(exon_mask,2)+1e-5)'+1e-5); % adapt number of reads starts according to profile
    seq_norm_weights = train_norm_sequence(CFG, seq_feat, seq_target);
    CFG.RR.seq_norm_weights = seq_norm_weights;
    norm_cov = norm_sequence(CFG, seq_feat); % predicted normalised coverage
    
    for n = 1:F*T,
      %fprintf('%i\r', n);
      exon_feat(:,n) = exon_feat(:,n) .* norm_cov; 
    end
    
    if 0
    cj = 0; ci = 0;
    norm_cov_sp = sparse(P_all, F*T);
    for g = 1:length(genes),
      for t = 1:length(genes(g).transcripts),
        norm_cov_sp(ci+[1:genes(g).exonic_len], cj+[1:F]) = repmat(norm_cov(ci+[1:genes(g).exonic_len]),1,F);
        cj = cj + F;
      end
      ci = ci + genes(g).exonic_len;
    end
    if exist(subs_idx, 'var')
      norm_cov_sp = norm_cov_sp(subs_idx, :);
    end
    exon_feat = exon_feat .* norm_cov_sp; % adapt exon features according to sequence normalisation
    end
    
    %exon_feat = exon_feat .* repmat(norm_cov, 1, F*T);
    
    %plot_sequence_weights(CFG.RR.seq_norm_weights, CFG.RR.order, 2*CFG.RR.half_win_size);
  end
  
  if CFG.norm_seqbias
    norm_weights = norm([weights_old, seq_norm_weights_old', reshape(profile_weights_old,1,F*N)] - [weights, seq_norm_weights', reshape(profile_weights,1,F*N)]);
  else
    norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N)] - [weights, reshape(profile_weights,1,F*N)]);
  end
  if fval_old(end)>=fval(end), sg = '-'; else sg = '+'; end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.3d\t%.3d\t%s\t%.1f\t%.1f\n', iter, fval(end), norm_weights, sg, 100/2*num_changed/sum([1:sum(pw_nnz)-1]), toc); end
  %if CFG.VERBOSE>1, fprintf(1, '%i\t%.3d\t%.3d\n', iter, fval(end), norm_weights); end
  if norm(fval_old-fval)<eps | norm_weights<eps | iter >= CFG.max_iter
    break;
  end
  %if iter>15, keyboard; end
  iter = iter + 1;
end
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

obj = fval(end);





