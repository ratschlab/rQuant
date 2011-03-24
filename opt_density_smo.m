function [profile_weights, obj, seq_weights] = opt_density_smo(CFG, genes)
% [profile_weights, obj, seq_weights] = opt_density_smo(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc.
%
% -- output --
% profile_weights: weights of profile functions 
% obj: objective evaluated with optimal parameters
% seq_weights: weights for sequence normalisation


DEBUG = 1;

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
coverage = sparse(P, 1);
exon_feat = sparse(P, F*T); % stores exon features from all transcripts in profile gene set
intron_count = zeros(I, 1);
intron_mask = zeros(I, T);
mask = true(P, 1); 
ci = 0; cj = 0; ct = 0; cn = 0;
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
  if DEBUG, fprintf('%i\r', g); end
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
    for t = 1:length(genes(g).transcripts),
      tmp_X = gen_sequence_features(CFG, genes(g), t);
      %[tmp_X tmp_Y tmp_Y_bg] = gen_sequence_features(CFG, genes(g), t);
      assert(size(tmp_X,1)==S);
      assert(size(tmp_X,2)==genes(g).exonic_len);
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
  ci = ci + genes(g).exonic_len;
  ct = ct + length(genes(g).transcripts);
  cn = cn + length(tmp_intron_count);
  clear tmp_introns tmp_intron_mask tmp_intron_count;
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
    seq_feat = seq_feat(subs_idx, :);
    %seq_target = seq_target(:, subs_idx);
    %seq_target_bg = seq_target_bg(:, subs_idx);
  end
  if CFG.VERBOSE>0, fprintf('subsampled from %i to %i positions\n', P_all, P); end
  clear P_old;
end

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
pw_nnz


% TODO
%   speed up
%%%%% optimisation %%%%%
eps = 1e-2;
C_w = [genes.transcript_length]';
%%% initialisation of variables
weights_old = zeros(1,T);
if CFG.norm_seqbias
  lambda = sum(sum(pw_nnz))/(sum(sum(pw_nnz))+S);
  seq_weights = (1-lambda)*ones(S, 1);
  seq_weights_old = zeros(S, 1);
  num_opt_steps = 1+sum([1:sum(sum(pw_nnz))-1])+sum([1:S-1]);
else
  lambda = 1;
  seq_weights = nan;
  num_opt_steps = 1+sum([1:sum(sum(pw_nnz))-1]);
end
profile_weights = zeros(F, N);
profile_weights(pw_nnz) = lambda*(F*N)/sum(sum(pw_nnz));
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
if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tDelta norm\n'); end
while 1
  tic
  weights_old = weights;
  fval_old = fval;
  profile_weights_old = profile_weights;
  if CFG.norm_seqbias
    seq_weights_old = seq_weights;
  end
  % pre-computations
  tmp_profiles = sparse(F*T,T);
  tmp_profiles(tp_idx) = profile_weights(:,tscp_len_bin);
  if CFG.norm_seqbias
    tmp_seq_weights = sparse(S*T,T);
    tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
  end
  
  if CFG.norm_seqbias
    seq_coeff = seq_feat*tmp_seq_weights;
  else
    seq_coeff = ones(P,T);
  end
  
  %%%%% A. optimise transcript weights
  exon_mask = exon_feat*tmp_profiles.*seq_coeff;
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  R_const = CFG.C_N*sum(sum((profile_weights(:,1:end-1)-profile_weights(:,2:end)).^2)) + CFG.C_F*sum(sum((profile_weights(1:end-1,:)-profile_weights(2:end,:)).^2));
  [weights, fval(1)] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, CFG.C_I, R_const, 1, weights, 'L1');
  assert(fval_old(end)-fval(1)>-1e-3);
  %corr(weights', [genes.expr_orig]')
  CFG.VERBOSE = tmp_VERBOSE;
    
  %%%%% B. optimise profile weights
  profile_weights = reshape(profile_weights, 1, F*N);
  R_const = abs(weights)*C_w;
  if size(intron_mask,1)>0, R_const = R_const + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
  num_changed = 0;
  pidx = find(pw_nnz(1:N*F-1));
  cnt = 1;
  ii = 0;
  for p = pidx,
    qidx = find(pw_nnz);
    qidx(qidx<=p) = [];
    ridx = randperm(length(qidx));
    qidx = qidx(ridx);
    for q = qidx,
      ii = ii + 1;
      if DEBUG, fprintf('%3.2f\r', 100*ii/sum([1:sum(pw_nnz)-1])); end
      theta1 = profile_weights(q);
      theta2 = profile_weights(p);
      d = theta1 + theta2;
      f1 = mod(q,F); if f1==0, f1 = F; end
      f2 = mod(p,F); if f2==0, f2 = F; end
      n1 = ceil(q/F);
      n2 = ceil(p/F);
      assert(f1~=f2 | n1~=n2);
      idx_w_n1 = find(tscp_len_bin==n1);  % all transcripts in length bin n1
      idx_w_n2 = find(tscp_len_bin==n2);  % all transcripts in length bin n2
      idx_w_th1 = f1+(idx_w_n1-1)*F;      % entries corresponding to theta1
      idx_w_th2 = f2+(idx_w_n2-1)*F;      % entries corresponding to theta2
      idx_wo_n1 = find(tscp_len_bin~=n1); % all transcripts not in n1
      idx_wo_n2 = find(tscp_len_bin~=n2); % all transcripts not in n2
      idx_wo_th1 = f1+(idx_wo_n1-1)*F;    % entries corresponding to theta_f1_notn1
      idx_wo_th2 = f2+(idx_wo_n2-1)*F;    % entries corresponding to theta_f2_notn2
      idx_wo_th12 = setdiff([1:F*T], [f1:F:F*T, f2:F:F*T]);
      %%% Rth2: residue for theta2
      Rth2 = exon_feat(:,idx_w_th2).*seq_coeff(:,idx_w_n2)*weights(idx_w_n2)' - exon_feat(:,idx_w_th1).*seq_coeff(:,idx_w_n1)*weights(idx_w_n1)';
      %%% R1: residue for theta1/theta2 independent variables
      R1 = exon_feat(:,idx_wo_th12)*tmp_profiles(idx_wo_th12,:).*seq_coeff*weights' + ...
           exon_feat(:,idx_w_th1).*seq_coeff(:,idx_w_n1)*weights(idx_w_n1)'*d - ...
           coverage;
      if f1==f2,
        idx_wo_th = intersect(idx_wo_th1, idx_wo_th2);
        idx_wo_n = intersect(idx_wo_n1, idx_wo_n2);
        R1 = R1 + exon_feat(:,idx_wo_th)*tmp_profiles(idx_wo_th,idx_wo_n).*seq_coeff(:,idx_wo_n)*weights(idx_wo_n)';
      else
        R1 = R1 + exon_feat(:,idx_wo_th1)*tmp_profiles(idx_wo_th1,idx_wo_n1).*seq_coeff(:,idx_wo_n1)*weights(idx_wo_n1)' + ...
             exon_feat(:,idx_wo_th2)*tmp_profiles(idx_wo_th2,idx_wo_n2).*seq_coeff(:,idx_wo_n2)*weights(idx_wo_n2)';
      end
      %%% R2: residue for coupling transcript length bins
      R2 = 0;
      if n1<N && ~(f1==f2 && n1+1==n2), R2 = R2 + (d-profile_weights(f1+n1*F))^2; end
      if n2<N
        if f1==f2 && n1==n2+1
          R2 = R2 + d^2;
        else
          R2 = R2 + profile_weights(f2+n2*F)^2;
        end
      end
      if n1>1 && ~(f1==f2 && n1==n2+1), R2 = R2 + (d-profile_weights(f1+(n1-2)*F))^2; end
      if n2>1
        if f1==f2 && n1+1==n2
          R2 = R2 + d^2;
        else
          R2 = R2 + profile_weights(f2+(n2-2)*F)^2;
        end
      end
      for f = 1:F,
        if f==f1 || f==f2, continue; end
        for n = 1:N-1,
          R2 = R2 + (profile_weights(f+(n-1)*F)-profile_weights(f+n*F))^2;
        end
      end
      for f = unique([f1 f2]),
        for n = 1:N-1,
          if f==f1 && (n==n1 || n==n1-1) || f==f2 && (n==n2 || n==n2-1), continue; end
          R2 = R2 + (profile_weights(f+(n-1)*F) - profile_weights(f+n*F))^2;
        end
      end
      %%% R3: residue for coupling supporting points
      R3 = 0;
      if f1<F && ~(n1==n2 && f1+1==f2), R3 = R3 + (d-profile_weights(f1+1+(n1-1)*F))^2; end
      if f2<F,
        if n1==n2 && f1==f2+1
          R3 = R3 + d^2;
        else
          R3 = R3 + profile_weights(f2+1+(n2-1)*F)^2;
        end
      end
      if f1>1 && ~(n1==n2 && f1==f2+1), R3 = R3 + (d-profile_weights(f1-1+(n1-1)*F))^2; end
      if f2>1,
        if n1==n2 && f1+1==f2
          R3 = R3 + d^2;
        else
          R3 = R3 + profile_weights(f2-1+(n2-1)*F)^2;
        end
      end
      for n = 1:N,
        if n==n1 || n==n2, continue; end
        for f = 1:F-1,
          R3 = R3 + (profile_weights(f+(n-1)*F) - profile_weights(f+1+(n-1)*F))^2;
        end
      end
      for n = unique([n1 n2]),
        for f = 1:F-1,
          if n==n1 && (f==f1 || f==f1-1) || n==n2 && (f==f2 || f==f2-1), continue; end
          R3 = R3 + (profile_weights(f+(n-1)*F) - profile_weights(f+1+(n-1)*F))^2;
        end
      end
      %%% S1: residue of quadratic term
      S1 = sum(Rth2.^2);
      if f1==f2
        if n1+1==n2 || n1==n2+1, S1 = S1 + CFG.C_N*6; else S1 = S1 + CFG.C_N*4; end
      else
        S1 = S1 + CFG.C_N*4;
      end
      if n1==n2
        if f1+1==f2 || f1==f2+1, S1 = S1 + CFG.C_F*6; else S1 = S1 + CFG.C_F*4; end
      else
        S1 = S1 + CFG.C_F*4;
      end
      if n1==1 || n1==N, S1 = S1 - CFG.C_N; end
      if n2==1 || n2==N, S1 = S1 - CFG.C_N; end
      if f1==1 || f1==F, S1 = S1 - CFG.C_F; end
      if f2==1 || f2==F, S1 = S1 - CFG.C_F; end
      assert(S1>0); % condition for minimum (2nd derivative > 0)
      %%% S2: residue of linear term
      S2 = sum(Rth2'*R1);
      % coupling constraints
      if n1<N && ~(f1==f2 && n1+1==n2), S2 = S2 + CFG.C_N*(profile_weights(f1+n1*F) - d); end       % theta_f1,n1+1
      if n1>1 && ~(f1==f2 && n1==n2+1), S2 = S2 + CFG.C_N*(profile_weights(f1+(n1-2)*F) -d); end    % theta_f1,n1-1
      if f1<F && ~(n1==n2 && f1+1==f2), S2 = S2 + CFG.C_F*(profile_weights(f1+1+(n1-1)*F) - d); end % theta_f1+1,n1
      if f1>1 && ~(n1==n2 && f1==f2+1), S2 = S2 + CFG.C_F*(profile_weights(f1-1+(n1-1)*F) - d); end % theta_f1-1,n1
      if n2<N, if f1==f2 && n1==n2+1, S2 = S2 - CFG.C_N*2*d; else S2 = S2 - CFG.C_N*profile_weights(f2+n2*F); end; end       % theta_f2,n2+1
      if n2>1, if f1==f2 && n1+1==n2, S2 = S2 - CFG.C_N*2*d; else S2 = S2 - CFG.C_N*profile_weights(f2+(n2-2)*F); end; end   % theta_f2,n2-1
      if f2<F, if n1==n2 && f1==f2+1, S2 = S2 - CFG.C_F*2*d; else S2 = S2 - CFG.C_F*profile_weights(f2+1+(n2-1)*F); end; end % theta_f2+1,n2
      if f2>1, if n1==n2 && f1+1==f2, S2 = S2 - CFG.C_F*2*d; else S2 = S2 - CFG.C_F*profile_weights(f2-1+(n2-1)*F); end; end % theta_f2-1,n2
      %%% S3: constant term
      S3 = sum(R1.^2) + CFG.C_N*R2 + CFG.C_F*R3 + R_const;
      %%% calculation and clipping of theta2 and theta1
      theta2_new = -S2/S1;
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
      end 
      if abs(theta2_new-theta2)>eps
        num_changed = num_changed + 1;
      end
      % update thetas and objective value
      profile_weights(q) = theta1_new;
      profile_weights(p) = theta2_new;
      cnt = cnt + 1;
      fval(cnt) = quad_fun(theta2_new, S1, 2*S2, S3);
      tmp_pw = reshape(profile_weights, F, N);
      tmp_profiles(tp_idx) = tmp_pw(:,tscp_len_bin);
      obj_alt = sum((exon_feat*tmp_profiles.*seq_coeff*weights'-coverage).^2) + R_const + CFG.C_N*sum(sum((tmp_pw(:,1:end-1)-tmp_pw(:,2:end)).^2)) + CFG.C_F*sum(sum((tmp_pw(1:end-1,:)-tmp_pw(2:end,:)).^2));
      assert(abs(fval(cnt)-obj_alt)<1e-3); % objective should be indentical to not-expanded objective
    end
  end
  assert(all(fval(1:cnt-1)-fval(2:cnt)>-1e-3)); % objective should decrease at every step
  profile_weights = reshape(profile_weights, F, N);
  %figure(iter); plot(profile_weights); ylim([0 10]);
  %plot(iter, fval(end), 'x');
  
  %%%%% C. determine sequence bias
  if CFG.norm_seqbias
    exon_mask = exon_feat*tmp_profiles;
    R_const = abs(weights)*C_w + CFG.C_N*sum(sum((tmp_pw(:,1:end-1)-tmp_pw(:,2:end)).^2)) + CFG.C_F*sum(sum((tmp_pw(1:end-1,:)-tmp_pw(2:end,:)).^2));
    if size(intron_mask,1)>0, R_const = R_const + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
    ii = 0;
    for k = 1:S-1,
      midx = k+1:S;
      %ridx = randperm(length(midx));
      %midx = midx(ridx);
      for m = midx,
        ii = ii + 1;
        if DEBUG, fprintf('%3.2f\r', 100*ii/sum([1:S-1])); end
        beta1 = seq_weights(m);
        beta2 = seq_weights(k);
        d = beta1 + beta2;
        %%% Rbe2: residue for beta2
        Rbe2 = exon_mask.*seq_feat(:,k:S:S*T)*weights' - exon_mask.*seq_feat(:,m:S:S*T)*weights';
        %%% R1: residue for beta1/beta2 independent variables
        idx_wo_be12 = setdiff([1:S*T], [m:S:S*T, k:S:S*T]);
        tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
        R1 = exon_mask.*(seq_feat(:,idx_wo_be12)*tmp_seq_weights(idx_wo_be12,:))*weights' + ...
             exon_mask.*seq_feat(:,m:S:S*T)*weights'*d - ...
             coverage;
        %%% S1: residue of quadratic term
        S1 = sum(Rbe2.^2);
        %%% S2: residue of linear term
        S2 = sum(Rbe2'*R1);
        %%% S3: constant term
        S3 = sum(R1.^2) + R_const;
        %%% calculation and clipping of beta2 and beta1
        beta2_new = -S2/S1;
        if beta2_new<0
          beta2_new = 0.0;
        end
        if beta2_new>d
          beta2_new = d;
        end
        beta1_new = d - beta2_new;
        % check if betas have been changed
        if abs(beta1_new-beta1)>eps
          num_changed = num_changed + 1;
        end 
        if abs(beta2_new-beta2)>eps
          num_changed = num_changed + 1;
        end
        % update betas and objective value
        seq_weights(m) = beta1_new;
        seq_weights(k) = beta2_new;
        cnt = cnt + 1;
        fval(cnt) = quad_fun(beta2_new, S1, 2*S2, S3);
        tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
        obj_alt = sum((exon_mask.*(seq_feat*tmp_seq_weights)*weights'-coverage).^2) + R_const;
        %%% CHECK THIS! %%%
        assert(abs(fval(cnt)-obj_alt)<1); % objective should be indentical to not-expanded objective
      end
    end
    assert(all(fval(1:cnt-1)-fval(2:cnt)>-1e-3)); % objective should decrease at every step
    %plot_sequence_weights(seq_weights, CFG.RR.order, 2*CFG.RR.half_win_size);
  end
  assert(cnt==num_opt_steps);
  
  if CFG.norm_seqbias
    norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N), seq_weights_old'] - [weights, reshape(profile_weights,1,F*N), seq_weights']);
  else
    norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N)] - [weights, reshape(profile_weights,1,F*N)]);
  end
  if fval_old(end)>=fval(end), sg = '-'; else sg = '+'; end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.3d\t%.3d\t%s\t%.1f\t%.1f\n', iter, fval(end), norm_weights, sg, 100/2*num_changed/num_opt_steps, toc); end
  
  %if CFG.VERBOSE>1, fprintf(1, '%i\t%.3d\t%.3d\n', iter, fval(end), norm_weights); end
  if norm(fval_old-fval)<eps || norm_weights<eps || iter >= CFG.max_iter
    break;
  end
  %if iter>15, keyboard; end
  iter = iter + 1;
end
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

obj = fval(end);





