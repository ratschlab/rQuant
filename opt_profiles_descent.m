function [profile_weights, obj, fval] = opt_profiles_descent(CFG, profile_weights, exon_feat, tmp_profiles, weights, coverage, seq_coeff, R_const, tscp_len_bin)
% [profile_weights, obj, fval] = opt_profiles_descent(CFG, profile_weights, exon_feat, tmp_profiles, weights, coverage, seq_coeff, R_const, tscp_len_bin)
%
% -- input --
% CFG: configuration struct
% profile_weights: weights of profile functions
% exon_feat: P x num_bins matrix of features for P exonic positions
% tmp_profiles: F*T x T matrix of profile weights for each transcript
% weights: weights of transcripts
% coverage: vector of observed exon coverage
% seq_coeff: vector of sequence correction
% R_const: constant residue
% tscp_len_bin: vector of length bins for each transcript
% 
%
% -- output --
% profile_weights: weights of profile functions
% obj: objective value at optimum
% fval: objective value at each step


F = CFG.num_plifs;                     % number of supporting points
N = size(CFG.transcript_len_ranges,1); % number of length bins
T = length(weights);
tp_idx = zeros(1,F*T);
for n = 1:F,
  tp_idx([1:F:length(tp_idx)]+n-1) = n+[0:(F*T+F):F*T*T];
end

max_iter = 1;

profile_weights = reshape(profile_weights, 1, F*N);
profile_weights_old = zeros(1, F*N);

pidx = find(profile_weights(1:N*F));

fval = 1e100*ones(1, length(pidx));
fval_old = zeros(1, length(pidx));

if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
iter = 1;
while 1
  cnt = 1;
  profile_weights_old = profile_weights;
  fval_old = fval;
  for p = pidx,
    f1 = mod(p,F); if f1==0, f1 = F; end
    n1 = ceil(p/F);
    idx_w_n = find(tscp_len_bin==n1);         % all transcripts in length bin n
    idx_w_th = f1+(idx_w_n-1)*F;              % entries corresponding to theta
    idx_wo_n = find(tscp_len_bin~=n1);        % all transcripts not in 1
    idx_w_th_wo_n = f1+(idx_wo_n-1)*F;        % entries corresponding to theta_f_notn
    idx_wo_th = setdiff([1:F*T], [f1:F:F*T]); % entries corresponding to theta_notf 
    %%% Rth: residue for theta
    if CFG.norm_seqbias
      Rth = exon_feat(:,idx_w_th).*seq_coeff(:,idx_w_n)*weights(idx_w_n)';
    else
      Rth = exon_feat(:,idx_w_th)*weights(idx_w_n)';
    end
    %%% R1: residue for theta independent variables
    if CFG.norm_seqbias
      R1 = exon_feat(:,idx_w_th_wo_n)*tmp_profiles(idx_w_th_wo_n,idx_wo_n).*seq_coeff(:,idx_wo_n)*weights(idx_wo_n)' + ...
           exon_feat(:,idx_wo_th)*tmp_profiles(idx_wo_th,:).*seq_coeff*weights' - ...
           coverage;
    else
      R1 = exon_feat(:,idx_w_th_wo_n)*tmp_profiles(idx_w_th_wo_n,idx_wo_n)*weights(idx_wo_n)' + ...
           exon_feat(:,idx_wo_th)*tmp_profiles(idx_wo_th,:)*weights' - ...
           coverage;
    end
    %%% R2: residue for coupling transcript length bins
    R2 = 0;
    if n1<N, R2 = R2 + profile_weights(f1+n1*F)^2; end
    if n1>1, R2 = R2 + profile_weights(f1+(n1-2)*F)^2; end
    for f = 1:F,
      for n = 1:N-1,
        if f~=f1 || (f==f1 && n~=n1 && n~=n1-1)
          R2 = R2 + (profile_weights(f+(n-1)*F)-profile_weights(f+n*F))^2;
        end
      end
    end
    %%% R3: residue for coupling supporting points
    R3 = 0;
    if f1<F, R3 = R3 + profile_weights(f1+1+(n1-1)*F)^2; end
    if f1>1, R3 = R3 + profile_weights(f1-1+(n1-1)*F)^2;; end
    for n = 1:N,
      for f = 1:F-1,
        if n~=n1 || (n==n1 && f~=f1 && f~=f1-1)
          R3 = R3 + (profile_weights(f+(n-1)*F)-profile_weights(f+1+(n-1)*F))^2;
        end
      end
    end
    %%% S1: residue of quadratic term
    S1 = sum(Rth.^2);
    % coupling constraints
    if n1<N, S1 = S1 + CFG.C_N; end
    if n1>1, S1 = S1 + CFG.C_N; end
    if f1<F, S1 = S1 + CFG.C_F; end
    if f1>1, S1 = S1 + CFG.C_F; end
    assert(S1>0); % condition for minimum (2nd derivative > 0)
    %%% S2: residue of linear term
    S2 = 2*sum(Rth'*R1);
    % coupling constraints
    if n1<N, S2 = S2 - 2*CFG.C_N*profile_weights(f1+n1*F); end
    if n1>1, S2 = S2 - 2*CFG.C_N*profile_weights(f1+(n1-2)*F); end
    if f1<F, S2 = S2 - 2*CFG.C_F*profile_weights(f1+1+(n1-1)*F); end
    if f1>1, S2 = S2 - 2*CFG.C_F*profile_weights(f1-1+(n1-1)*F); end
    %%% S3: constant term
    S3 = sum(R1.^2) + CFG.C_N*R2 + CFG.C_F*R3 + R_const;
    %%% calculation and clipping of theta
    th_new = -0.5*S2/S1;
    if th_new < 0
      profile_weights(p) = 0.0;
    else
      profile_weights(p) = th_new;
    end
    fval(cnt) = quad_fun(profile_weights(p), S1, S2, S3);
    tmp_pw = reshape(profile_weights, F, N);
    tmp_profiles(tp_idx) = tmp_pw(:,tscp_len_bin);
    if CFG.norm_seqbias
      obj_alt = sum((exon_feat*tmp_profiles.*seq_coeff*weights'-coverage).^2) + R_const + CFG.C_N*sum(sum((tmp_pw(:,1:end-1)-tmp_pw(:,2:end)).^2)) + CFG.C_F*sum(sum((tmp_pw(1:end-1,:)-tmp_pw(2:end,:)).^2));
    else
      obj_alt = sum((exon_feat*tmp_profiles*weights'-coverage).^2) + R_const + CFG.C_N*sum(sum((tmp_pw(:,1:end-1)-tmp_pw(:,2:end)).^2)) + CFG.C_F*sum(sum((tmp_pw(1:end-1,:)-tmp_pw(2:end,:)).^2));
    end
    assert(abs(fval(cnt)-obj_alt)<1e-3); % objective should be indentical to not-expanded objective
    cnt = cnt + 1;
  end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\n', iter, fval(end), norm(profile_weights_old-profile_weights)); end
  if norm(fval_old-fval)<1e-5 || norm(profile_weights_old-profile_weights)<1e-5 || iter>=max_iter,
    break;
  end
  iter = iter + 1;
end
assert(all(fval(1:end-1)-fval(2:end)>-1e-3)); % objective should decrease at every step
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

profile_weights = reshape(profile_weights, F, N);
obj = fval(end);