function [profile_weights, obj] = opt_profiles_smo(CFG, genes)

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


%%%%% optimisation
C_w = [genes.transcript_length]';
% initialisation of variables
weights_old = zeros(1,T);
profile_weights = ones(F, size(CFG.transcript_len_ranges,1));
profile_weights_old = zeros(F, size(CFG.transcript_len_ranges,1));
fval = 1e100; %1e100*ones(1,T+F*N);
fval_old = 0;
iter = 1;
if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
%figure(1); hold on;
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
  [weights, fval] = opt_transcripts_L2(CFG, coverage, exon_mask, [], [], C_w, 1);
  CFG.VERBOSE = tmp_VERBOSE;
    
  % B. optimise profile weights
  profile_weights = reshape(profile_weights, 1, F*N);
  num_changed = 0; examine_all = true;
  % TODO: clever choice of theta1 and theta2
    for p = 1:2:N*F-1,
    theta1 = profile_weights(p);
    theta2 = profile_weights(p+1);
    d = theta1 + theta2;
    f1 = mod(p,F); if f1==0, f1 = F; end
    f2 = mod(p+1,F); if f2==0, f2 = F; end
    n1 = ceil(p/F);
    n2 = ceil((p+1)/F);
    idx_w_n1 = find(tscp_len_bin==n1);  % all transcripts in length bin n1
    idx_w_n2 = find(tscp_len_bin==n2);  % all transcripts in length bin n2
    idx_w_th1 = f1+(idx_w_n1-1)*F;      % entries corresponding to theta1
    idx_w_th2 = f2+(idx_w_n2-1)*F;      % entries corresponding to theta2
    idx_wo_n1 = find(tscp_len_bin~=n1); % all transcripts not in n1
    idx_wo_n2 = find(tscp_len_bin~=n2); % all transcripts not in n2
    idx_wo_th1 = f1+(idx_wo_n1-1)*F;    % entries corresponding to theta_f1_notn1
    idx_wo_th2 = f2+(idx_wo_n2-1)*F;    % entries corresponding to theta_f2_notn2
    idx_wo_f = setdiff([1:F*T], [f1:F:F*T, f2:F:F*T]);
    % residue for theta2
    Rth2 = exon_feat(:,idx_w_th2)*tmp_profiles(idx_w_th2,idx_w_n2)*weights(idx_w_n2)' - ...
           exon_feat(:,idx_w_th1)*tmp_profiles(idx_w_th1,idx_w_n1)*weights(idx_w_n1)';
    % residue for theta1/theta2 independent variables
    Rp = exon_feat(:,idx_wo_th1)*tmp_profiles(idx_wo_th1,idx_wo_n1)*weights(idx_wo_n1)' + ...
         exon_feat(:,idx_wo_th2)*tmp_profiles(idx_wo_th2,idx_wo_n2)*weights(idx_wo_n2)' + ...
         exon_feat(:,idx_wo_f)*tmp_profiles(idx_wo_f,:)*weights' + ...
         exon_feat(:,idx_w_th1)*tmp_profiles(idx_w_th1,idx_w_n1)*weights(idx_w_n1)'*d - ...
         coverage;
    S1 = 2*sum(Rth2.^2) + 8; % quadratic term
    assert(S1>0); % condition for minimum (2nd derivative > 0)
    S2 = 4*d - sum(Rth2'*Rp); % linear term
    S3 = 0; % constant term -- to implement
    % coupling constraints
    if n1<N, S2 = S2 - 2*profile_weights(f1+n1*F); end % theta_f1,n1+1
    if n2<N, S2 = S2 + 2*profile_weights(f2+n2*F); end % theta_f2,n2+1
    if f1<F, S2 = S2 - 2*profile_weights(f1+1+(n1-1)*F); end % theta_f1+1,n1
    if f2<F, S2 = S2 + 2*profile_weights(f2+1+(n2-1)*F); end % theta_f2+1,n2
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
    eps = 1e^-3;
    if abs(theta1_new-theta1)>eps | abs(theta2_new-theta2)>eps
      num_changed = num_changed + 1;
    end
    % update thetas and objective value
    profile_weights(p) = theta1_new;
    profile_weights(p+1) = theta2_new;
    fval = quad_fun(theta2_new, S1/2, -S2, S3);
    tmp_pw = reshape(profile_weights, F, N);
    tmp_profiles(idx) = tmp_pw(:,tscp_len_bin);
  end
  profile_weights = reshape(profile_weights, F, N);
  %figure(iter); plot(profile_weights); ylim([0 10]);
  %plot(iter, fval(end), 'x');
  %if iter>20, keyboard; end
  norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N)] - [weights, reshape(profile_weights,1,F*N)]);
  if fval_old(end)>=fval(end), sg = '-'; else sg = '+'; end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\t%s\t%.1f\n', iter, fval(end), norm_weights, sg, 2*100*num_changed/(F*N)); end
  if norm(fval_old-fval)<1e-5 | norm_weights<1e-5
    break;
  end
  iter = iter + 1;
  %keyboard
end
obj = fval(end);
keyboard





