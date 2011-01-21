function [profile_weights, obj] = opt_profiles_smo(CFG, genes)
% [profile_weights, obj] = opt_profiles_smo(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc. 
%
% -- output --
% profile_weights: weights of profile functions 
% obj: objective evaluated with optimal parameters


T = length([genes.transcripts]);       % number of transcripts
P = sum([genes.exonic_len]);           % number of positions
F = CFG.num_plifs;                     % number of supporting points
N = size(CFG.transcript_len_ranges,1); % number of length bins


%%%%% pre-processing
exon_feat = sparse(P, F*T); % stores exon features from all transcripts in profile gene set
coverage = sparse(P, 1);
ci = 0; cj = 0; ct = 0;
weights = zeros(1,T);
for g = 1:length(genes),
  for t = 1:length(genes(g).transcripts),
    tmp_feat = gen_exon_features(genes(g), t, F, CFG.max_side_len, 1);
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
    exon_feat(ci+[1:genes(g).exonic_len], cj+[1:F]) = tmp_feat;
    cj = cj + F;
  end
  try
    [tmp_coverage excluded_reads reads_ok] = get_coverage_per_read(CFG, genes(g), 1);
    coverage(ci+[1:genes(g).exonic_len],1) = tmp_coverage;
  catch
    reads_ok = 0;
  end
  Tg = length(genes(g).transcripts);
  % initialisation of transcript weights to proportionate mean coverage
  weights(ct+[1:Tg]) = full(mean(coverage(ci+[1:genes(g).exonic_len]))/Tg*ones(1,Tg));
  ci = ci + genes(g).exonic_len;
  ct = ct + length(genes(g).transcripts);
  %if ~reads_ok, continue; end  
end
tscp_len_bin = [genes.transcript_len_bin];
% find thetas that do not need to be optimised (located in the body of the profile function)
lmt = linspace(0, sqrt(CFG.max_side_len), (F/2)+1).^2;
lmt(end) = inf;
pw_nnz = false(F, N);
max_tl = max([genes.transcript_length]);
if max_tl<CFG.transcript_len_ranges(N,1)
  CFG.transcript_len_ranges(N,2) = CFG.transcript_len_ranges(N,1) + 1;
else
  CFG.transcript_len_ranges(N,2) = max_tl;
end
for n = 1:N,
  fidx = find(CFG.transcript_len_ranges(n,2)/2>lmt, 1, 'last');
  pw_nnz([1:fidx, F-fidx+1:F],n) = true;
end
pw_nnz = reshape(pw_nnz, 1, F*N);


%%%%% optimisation
max_iter = 200;
C_w = [genes.transcript_length]';
% initialisation of variables
weights_old = zeros(1,T);
profile_weights = ones(F, N);
profile_weights_old = zeros(F, N);
fval = 1e100; %1e100*ones(1,T+F*N);
fval_old = 0;
iter = 1;
if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
while 1
  weights_old = weights;
  fval_old = fval;
  profile_weights_old = profile_weights;
  % pre-computations
  tmp_profiles = sparse(zeros(F*T,T));
  idx = zeros(1,F*T);
  for n = 1:F,
    idx([1:F:length(idx)]+n-1) = n+[0:(F*T+F):F*T*T];
  end
  tmp_profiles(idx) = profile_weights(:,tscp_len_bin);
  
  % A. optimise transcript weights
  exon_mask = exon_feat*tmp_profiles;
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  [weights, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, [], [], C_w, 1, weights, 'L1'); % intron model missing -- to implement
  %corr(weights', [genes.expr_orig]')
  CFG.VERBOSE = tmp_VERBOSE;
    
  % B. optimise profile weights
  profile_weights = reshape(profile_weights, 1, F*N);
  num_changed = 0; examine_all = true;
  changed = zeros(1, F*N);
  % TODO: clever choice of theta1 and theta2
  pidx = find(pw_nnz(1:N*F-1));
  ii = 0;
  for p = pidx,
    qidx = find(pw_nnz); 
    qidx(qidx<=p) = [];
    ridx = randperm(length(qidx));
    qidx = qidx(ridx);
    for q = qidx%p+1:N*F,
      ii = ii + 1;
      fprintf('%3.2f\r', 100*ii/sum([1:sum(pw_nnz)-1]));
      theta1 = profile_weights(p);
      theta2 = profile_weights(q);
      d = theta1 + theta2;
      f1 = mod(p,F); if f1==0, f1 = F; end
      f2 = mod(q,F); if f2==0, f2 = F; end
      n1 = ceil(p/F);
      n2 = ceil((q)/F);
      if n1~=n2, continue; end
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
      eps = 1e-5;
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
      tmp_profiles(idx) = tmp_pw(:,tscp_len_bin);
    end
  end
  profile_weights = reshape(profile_weights, F, N);
  %figure(iter); plot(profile_weights); ylim([0 10]);
  %plot(iter, fval(end), 'x');
  norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N)] - [weights, reshape(profile_weights,1,F*N)]);
  if fval_old(end)>=fval(end), sg = '-'; else sg = '+'; end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\t%s\t%.1f\n', iter, fval(end), norm_weights, sg, 100/2*num_changed/sum([1:sum(pw_nnz)-1])); end
  if norm(fval_old-fval)<1e-5 | norm_weights<1e-5 | iter >= max_iter
    break;
  end
  %if iter>15, keyboard; end
  iter = iter + 1;
  %keyboard
end
obj = fval(end);





