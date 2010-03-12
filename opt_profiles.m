function [weights, dists, genes] = opt_profiles(CFG, genes)
% [weights, dists, genes] = opt_profiles(CFG, genes)
%
% -- input --
% CFG: configuration struct
% genes: struct defining genes with start, stops, exons etc. 
%
% -- output --
% weights: weights of profile functions 
% dists: distances to closest intron
% genes: struct defining genes with respective expression and
% transcript length bin added


%%%%% define optimisation problem %%%%%
% min_{x=[w,d,xis]}  sum_{e=1}^{E} sum_{n=1}^{N} sum_{f=1}^{F} w_{e,n,f}
%                  + sum_{l=1}^{D} sum_{r=1}^{D} w_{l,r}
%                  + 0.5 sum xi^2
%                  + 0.5 C_tau sum tau^2
%                  + 0.5 C_kappa sum kappa^2
%                  + 0.5 C_kappa sum theta^2
%
% s.t.
% for all p=1..P:   -sum_{t=1}^{T} sum_{e=1}^{E} sum_{n=1}^{N}
%                    sum_{f=1}^{F} [[delta(len_{t})=n v
%                                    epsilon(w_{t}^{~}=e v
%                                    phi(t,p)=f)]] w_{e,n,f}
%                   + sum_{l=1}^{D} sum_{r=1}^{D} w_{l,r}
%                    - xi_{p}                                  = -ec_{p}
%                    sum_{e,n,f} 1/(E*N*F) w_{e,n,f}           = 1
% for all e=1..E,    
%         n=1..N,    w_{e,n,f}-w_{e,n,f+1} <= tau_{e,n,f}
%         f=1..F-1: -w_{e,n,f}+w_{e,n,f+1} <= tau_{e,n,f}
%                   
% for all e=1..E,
%         n=1..N-1   w_{e,n,f}-w_{e,n+1,f} <= kappa^{1}_{e,n,f}
%         f=1..F:   -w_{e,n,f}+w_{e,n+1,f} <= kappa^{1}_{e,n,f}
%                   
% for all e=1..E-1,
%         n=1..N,    w_{e,n,f}-w_{e+1,n,f} <= kappa^{2}_{e,n,f}
%         f=1..F:   -w_{e,n,f}+w_{e+1,n,f} <= kappa^{2}_{e,n,f}    
%                   
% for all l=1..D-1,  w_{l,r}-w_{l+1,r} <= theta^{1}_{l,r}
%         r=1..D:   -w_{l,r}+w_{l+1,r} <= theta^{1}_{l,r} 
%                   
% for all l=1..D,    w_{l,r}-w_{l,r+1} <= theta^{2}_{l,r}
%         r=1..D-1: -w_{l,r}+w_{l,r+1} <= theta^{2}_{l,r} 
% 
% with bounds:
%    0 <= w_{e,n,f}         <= INF  for all e=1..E, n=1..N, f=1..F
%    0 <= tau_{e,n,f}       <= INF  for all e=1..E, n=1..N, f=1..F-1
%    0 <= kappa^{1}_{e,n,f} <= INF  for all e=1..E, n=1..N-1, f=1..F
%    0 <= kappa^{2}_{e,n,f} <= INF  for all e=1..E-1, n=1..N, f=1..F
%    0 <= theta^{1}_{l,r}   <= INF  for all l=1..D-1, r=1..D
%    0 <= theta^{2}_{l,r}   <= INF  for all l=1..D, r=1..D-1
% -INF <= xi_{p}            <= INF  for all p in P
%
% P: number of positions
% F: number of supporting points of profile plifs
% N: number of transcript length bins
% E: number of expression bins
% D: number of supporting points of plifs for intron distances


INF = 1e20;


%%%%% preprocessing %%%%%
P = sum([genes.exonic_len]);
TR.mask = true(P,1);
cnt = 1; nnz_feat = 0;
for g = 1:length(genes),
  gene = genes(g);
  exon_mask = zeros(1, gene.exonic_len);
  exon_mask_expr = zeros(1, gene.exonic_len);
  dist = zeros(length(genes(g).transcripts), gene.exonic_len, CFG.num_intron_plifs, CFG.num_intron_plifs);
  for t = 1:length(genes(g).transcripts),
    [tmp_feat, genes(g).transcript_length(t)] = gen_exon_features(gene, t, CFG.num_plifs, CFG.max_side_len, 1);
    exon_mask = exon_mask + sum(tmp_feat>0,2)';
    exon_mask_expr = exon_mask_expr + sum(tmp_feat>0,2)'*(gene.transcript_weights(t)>0);
    dist(t,:,:,:) = gen_intron_features(gene, t, CFG.num_intron_plifs, CFG.read_len, 1)*gene.transcript_weights(t);
  end
  if CFG.subsample,
    subsample_mask = (rand(sum(exon_mask>0),1)<CFG.subsample_frac);
  else
    subsample_mask = true(sum(exon_mask>0),1);
  end
  TR.mask(cnt:cnt+sum(exon_mask>0)-1) = TR.mask(cnt:cnt+sum(exon_mask>0)-1) & subsample_mask;
  tmp = exon_mask_expr(exon_mask>0);
  nnz_feat = nnz_feat + sum(tmp(TR.mask(cnt:cnt+sum(exon_mask>0)-1)));
  dist = squeeze(sum(dist,1));
  tmp = dist(exon_mask>0,:,:);
  nnz_feat = nnz_feat + sum(sum(sum(tmp(TR.mask(cnt:cnt+sum(exon_mask>0)-1),:,:)>0,3),2),1);
  cnt = cnt + sum(exon_mask>0);
end
assert((cnt-1)==P); assert(P==size(TR.mask,1)); 

if CFG.subsample,
  fprintf('subsampled from %i to %i positions\n', P, sum(TR.mask));
end
P = sum(TR.mask);


%%%%% preparation of optimisation problem %%%%%
F = CFG.num_plifs;                     % number of supporting points of profile plifs
N = size(CFG.transcript_len_ranges,1); % number of transcript length bins
E = size(CFG.expr_ranges,1);           % number of expression bins
D = CFG.num_intron_plifs;              % number of supporting points of plifs for intron distances

tau_len = E*N*(F-1);    % taus couple adjacent supporting points 
kappa1_len = E*(N-1)*F; % kappa1s couple adjacent transcript length bins
kappa2_len = (E-1)*N*F; % kappa2s couple adjacent expression bins
kappa_len = kappa1_len + kappa2_len;
theta1_len = D*(D-1);
theta2_len = D*(D-1);
theta_len = theta1_len + theta2_len; % thetas

b = zeros(P + 1 + 2*tau_len + 2*kappa_len, 1);
Ac = 0; % matrix counter 
Ai = zeros(1, P + nnz_feat + E*N*F + 6*tau_len + 6*kappa_len + 6*theta_len); % row indices
Aj = zeros(1, P + nnz_feat + E*N*F + 6*tau_len + 6*kappa_len + 6*theta_len); % column indices
Av = zeros(1, P + nnz_feat + E*N*F + 6*tau_len + 6*kappa_len + 6*theta_len); % matrix values

p_offset = 0;
mask_offset = 0;
cnt = 0;
for g = 1:length(genes),
  g
  gene = genes(g);
  exon_mask = zeros(1, gene.exonic_len);
  profiles = zeros(length(genes(g).transcripts), gene.exonic_len, CFG.num_plifs);
  dist = zeros(length(genes(g).transcripts), gene.exonic_len, CFG.num_intron_plifs, CFG.num_intron_plifs);
  for t = 1:length(genes(g).transcripts),
    [tmp_feat, genes(g).transcript_length(t)] = gen_exon_features(gene, t, CFG.num_plifs, CFG.max_side_len, 1);
    profiles(t,:,:) = tmp_feat*gene.transcript_weights(t);
    exon_mask = exon_mask + sum(tmp_feat>0,2)';
    genes(g).expr_bin(t) = find(CFG.expr_ranges(:,1) <= gene.transcript_weights(t) & ...
                                CFG.expr_ranges(:,2) >= gene.transcript_weights(t));
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
    dist(t,:,:,:) = gen_intron_features(gene, t, CFG.num_intron_plifs, CFG.read_len, 1)*gene.transcript_weights(t);
  end
  assert(max(exon_mask)<=length(genes(g).transcripts));
  exon_mask_idx = find(exon_mask>0);
  profiles = profiles(:,exon_mask_idx,:);
  profiles_len = size(profiles,2);
  profiles = profiles(:,TR.mask([1:profiles_len]+mask_offset)>0,:);
  dist = squeeze(sum(dist(:,exon_mask_idx,:,:),1));
  dist = dist(TR.mask([1:profiles_len]+mask_offset)>0,:,:);
  assert(all(sum(sum(dist,3),2)>=0));
  %%%%% exon count constraints %%%%%
  for p = 1:size(profiles,2),
    % taus and kappas
    for t = 1:length(genes(g).transcripts),
      sp_idx = find(profiles(t,p,:)~=0);
      assert(length(sp_idx)<=1);
      cnt = cnt + length(sp_idx);
      Ai(Ac+1:Ac+length(sp_idx)) = p_offset + p;
      Aj(Ac+1:Ac+length(sp_idx)) = sp_idx + (gene.expr_bin(t)-1)*N*F+(genes(g).transcript_len_bin(t)-1)*F;
      Av(Ac+1:Ac+length(sp_idx)) = -profiles(t, p, sp_idx); 
      Ac            = Ac + length(sp_idx);
    end
    % thetas
    tmp_dist = squeeze(dist(p,:,:));
    sp_idx = find(tmp_dist~=0);
    assert(~isempty(sp_idx) | (any(genes(g).transcript_weights==0) & isempty(sp_idx)));
    cnt = cnt + length(sp_idx);
    Ai(Ac+1:Ac+length(sp_idx)) = p_offset + p;
    Aj(Ac+1:Ac+length(sp_idx)) = E*N*F + sp_idx;
    Av(Ac+1:Ac+length(sp_idx)) = tmp_dist(sp_idx);
    Ac            = Ac + length(sp_idx);
    % xis
    Ai(Ac+1)      = p_offset + p;
    Aj(Ac+1)      = E*N*F + D*D + p_offset + p;
    Av(Ac+1)      = -1;
    Ac            = Ac + 1;
  end
  try
    [coverage excluded_reads pair_ok] = get_coverage_per_read(CFG, gene, 1);
    %if CFG.norm_seqbias
    %  coverage = norm_sequence(CFG, gene, coverage);
    %end
    coverage = sum(coverage,2);
  catch
    pair_ok = 0;
  end;
  if ~pair_ok, continue; end
  b([1:size(profiles,2)]+p_offset,1) = -coverage(exon_mask_idx(TR.mask([1:profiles_len]+mask_offset)>0));
  mask_offset = mask_offset + profiles_len;
  p_offset = p_offset + size(profiles,2);
end
assert(p_offset==P);
assert(nnz_feat==cnt);

clear coverage tmp_feat profile;


%%%%% constraint: sum of weights = 1 %%%%%
offset = P;
% sum_f w_f/F = 1 
offset = offset + 1;
Ai(Ac+[1:E*N*F]) = offset;
Aj(Ac+[1:E*N*F]) = 1:E*N*F;
Av(Ac+[1:E*N*F]) = 1/(E*N*F); 
b(offset, 1) = 1;
Ac = Ac + E*N*F;


%%%%% tau contraints: coupling adjacent supporting points %%%%%
%%%%% |w_{e,n,f}-w{e,n,f+1}|=tau_{e,n,f}                  %%%%%
for e = 1:E,
  for n = 1:N,
    for f = 1:F-1,
      % w_f - w_{f+1} = tau_f
      % w_f
      rc = offset + 2*f-1 + (e-1)*N*2*(F-1) + (n-1)*2*(F-1); % index to row number
      Ai(Ac+1) = rc;
      Aj(Ac+1) = f + (e-1)*N*F + (n-1)*F;
      Av(Ac+1) = 1; 
      % -w_{f+1}
      Ai(Ac+2) = rc;
      Aj(Ac+2) = f + 1 + (e-1)*N*F + (n-1)*F;
      Av(Ac+2) = -1;
      % -tau_f
      Ai(Ac+3) = rc;
      Aj(Ac+3) = E*N*F + D*D + P + f + (e-1)*N*(F-1) + (n-1)*(F-1);
      Av(Ac+3) = -1;
      b(rc, 1) = 0;
      Ac       = Ac + 3;
      % -w_f + w_{f+1} = tau_f
      % -w_f
      Ai(Ac+1) = rc + 1;
      Aj(Ac+1) = f + (e-1)*N*F + (n-1)*F;
      Av(Ac+1) = -1; 
      % w_{f+1}
      Ai(Ac+2) = rc + 1;
      Aj(Ac+2) = f + 1 + (e-1)*N*F + (n-1)*F;
      Av(Ac+2) = 1;
      % -tau_f
      Ai(Ac+3) = rc + 1;
      Aj(Ac+3) = E*N*F + D*D + P + f + (e-1)*N*(F-1) + (n-1)*(F-1);
      Av(Ac+3) = -1;
      b(rc+1, 1) = 0;
      Ac       = Ac + 3;
    end
  end 
end


%%%%% kappa contraints %%%%%
for f = 1:F,
  %%%%% kappa1 contraints: coupling adjacent transcript length bins %%%%%
  %%%%% |w_{e,n,f}-w{e,n+1,f}|=kappa^{1}_{e,n,f}                    %%%%%
  for e = 1:E,
    for n = 1:N-1,
      rc =  offset + 2*tau_len + 2*f-1 + (e-1)*(N-1)*2*F + (n-1)*2*F; % index to row number
      % w_n
      Ai(Ac+1) = rc;
      Aj(Ac+1) = f + (e-1)*N*F + (n-1)*F;
      Av(Ac+1) = 1; 
      % -w_{n+1}
      Ai(Ac+2) = rc;
      Aj(Ac+2) = f + (e-1)*N*F + (n-1+1)*F;
      Av(Ac+2) = -1;
      % -kappa1_n
      Ai(Ac+3) = rc;
      Aj(Ac+3) = E*N*F + D*D + P + tau_len + f + (e-1)*(N-1)*F + (n-1)*F;
      Av(Ac+3) = -1;
      b(rc, 1) = 0;
      Ac       = Ac + 3;
      % -w_n
      Ai(Ac+1) = rc + 1;
      Aj(Ac+1) = f + (e-1)*N*F + (n-1)*F;
      Av(Ac+1) = -1; 
      % w_{n+1}
      Ai(Ac+2) = rc + 1;
      Aj(Ac+2) = f + (e-1)*N*F + (n-1+1)*F;
      Av(Ac+2) = 1;
      % -kappa1_n
      Ai(Ac+3) = rc + 1;
      Aj(Ac+3) = E*N*F + D*D + P + tau_len + f + (e-1)*(N-1)*F + (n-1)*F;
      Av(Ac+3) = -1;
      b(rc+1, 1) = 0;
      Ac       = Ac + 3;
    end 
  end
  %%%%% kappa2 contraints: coupling adjacent expression bins %%%%%
  %%%%% |w_{e,n,f}-w{e+1,n,f}|=kappa^{2}_{e,n,f}             %%%%%
  for n = 1:N,
    for e = 1:E-1,
      rc =  offset + 2*tau_len + 2*kappa1_len + 2*f-1 + (e-1)*N*2*F + (n-1)*2*F; % index to row number
      % w_e
      Ai(Ac+1) = rc;
      Aj(Ac+1) = f + (e-1)*N*F + (n-1)*F;
      Av(Ac+1) = 1; 
      % -w_{e+1}
      Ai(Ac+2) = rc;
      Aj(Ac+2) = f + (e-1+1)*N*F + (n-1)*F;
      Av(Ac+2) = -1;
      % -kappa2_e
      Ai(Ac+3) = rc;
      Aj(Ac+3) = E*N*F + D*D + P + tau_len + kappa1_len + f + (e-1)*N*F + (n-1)*F;
      Av(Ac+3) = -1;
      b(rc, 1) = 0;
      Ac       = Ac + 3;
      % -w_e
      Ai(Ac+1) = rc + 1;
      Aj(Ac+1) = f + (e-1)*N*F + (n-1)*F;
      Av(Ac+1) = -1; 
      % w_{e+1}
      Ai(Ac+2) = rc + 1;
      Aj(Ac+2) = f + (e-1+1)*N*F + (n-1)*F;
      Av(Ac+2) = 1;
      % -kappa2_e
      Ai(Ac+3) = rc + 1;
      Aj(Ac+3) = E*N*F + D*D + P + tau_len + kappa1_len + f + (e-1)*N*F + (n-1)*F;
      Av(Ac+3) = -1;
      b(rc+1, 1) = 0;
      Ac       = Ac + 3;
    end 
  end
end


%%%%% theta contraints %%%%%
%%%%% theta1 contraints: coupling adjacent supporting points for 3' introns %%%%%
%%%%% |w_{l,r}-w{l+1,r}|=theta^{1}_{l,r}                                    %%%%%
for r = 1:D,
  for l = 1:D-1,
    % w_l - w_{l+1} = theta1_l
    % w_l
    rc = offset + 2*tau_len + 2*kappa_len + 2*l-1 + (r-1)*2*(D-1); % index to row number
    Ai(Ac+1) = rc;
    Aj(Ac+1) = E*N*F + l + (r-1)*D;
    Av(Ac+1) = 1; 
    % -w_{l+1}
    Ai(Ac+2) = rc;
    Aj(Ac+2) = E*N*F + l + 1 + (r-1)*D;
    Av(Ac+2) = -1;
    % -theta1_l
    Ai(Ac+3) = rc;
    Aj(Ac+3) = E*N*F + D*D + P + tau_len + kappa_len + l + (r-1)*(D-1);
    Av(Ac+3) = -1;
    b(rc, 1) = 0;
    Ac       = Ac + 3;
    % -w_l + w_{l+1} = theta1_l
    % -w_l
    Ai(Ac+1) = rc + 1;
    Aj(Ac+1) = E*N*F + l + (r-1)*D;
    Av(Ac+1) = -1; 
    % w_{l+1}
    Ai(Ac+2) = rc + 1;
    Aj(Ac+2) = E*N*F + l + 1 + (r-1)*D;
    Av(Ac+2) = 1;
    % -theta1_l
    Ai(Ac+3) = rc + 1;
    Aj(Ac+3) = E*N*F + D*D + P + tau_len + kappa_len + l + (r-1)*(D-1);
    Av(Ac+3) = -1;
    b(rc+1, 1) = 0;
    Ac       = Ac + 3;
  end 
end
%%%%% theta2 contraints: coupling adjacent supporting points for 5' introns %%%%%
%%%%% |w_{l,r}-w{l,r+1}|=theta^{2}_{l,r}                                    %%%%%
for l = 1:D,
  for r = 1:D-1,
    % w_r - w_{r+1} = theta2_r
    % w_r
    rc = offset + 2*tau_len + 2*kappa_len + 2*theta1_len + 2*r-1 + (l-1)*2*(D-1); % index to row number
    Ai(Ac+1) = rc;
    Aj(Ac+1) = E*N*F + l + (r-1)*D;
    Av(Ac+1) = 1; 
    % -w_{r+1}
    Ai(Ac+2) = rc;
    Aj(Ac+2) = E*N*F + l + (r-1+1)*D;
    Av(Ac+2) = -1;
    % -theta2_r
    Ai(Ac+3) = rc;
    Aj(Ac+3) = E*N*F + D*D + P + tau_len + kappa_len + theta1_len + r + (l-1)*(D-1);
    Av(Ac+3) = -1;
    b(rc, 1) = 0;
    Ac       = Ac + 3;
    % -w_r + w_{r+1} = theta2_r
    % -w_r
    Ai(Ac+1) = rc + 1;
    Aj(Ac+1) = E*N*F + l + (r-1)*D;
    Av(Ac+1) = -1; 
    % w_{r+1}
    Ai(Ac+2) = rc + 1;
    Aj(Ac+2) = E*N*F + l + 1 + (r-1)*D;
    Av(Ac+2) = 1;
    % -theta2_r
    Ai(Ac+3) = rc + 1;
    Aj(Ac+3) = E*N*F + D*D + P + tau_len + kappa_len + theta1_len + r + (l-1)*(D-1);
    Av(Ac+3) = -1;
    b(rc+1, 1) = 0;
    Ac       = Ac + 3;
  end 
end


assert(all(Ai>0 & max(Ai)==offset+2*(tau_len+kappa_len+theta_len)) & all(Aj>0 & max(Aj)==E*N*F+D*D+P+tau_len+kappa_len+theta_len));
% (P+1+2*(tau_len+kappa_len+theta_len))x(E*N*F+D*D+P+tau_len+kappa_len+theta_len) matrix
A = sparse(Ai, Aj, Av, offset + 2*tau_len + 2*kappa_len + 2*theta_len , E*N*F + D*D + P + tau_len + kappa_len + theta_len);
clear Ai Aj Av feat;


%%%%% bounds %%%%%
LB = [zeros(E*N*F+D*D,1); -INF*ones(P,1); zeros(tau_len+kappa_len+theta_len, 1)]; % lower bounds for xopt
UB = INF*ones(E*N*F+D*D+P+tau_len+kappa_len+theta_len,1); % upper bounds for xopt

if CFG.num_intron_plifs>2,
  LB(E*N*F + D*D) = 0;
  UB(E*N*F + D*D) = 0;
else % no intron distance constraints
  LB(E*N*F+1:E*N*F + D*D) = 0;
  UB(E*N*F+1:E*N*F + D*D) = 0;
end


%%%%% linear term of the objective function %%%%%
obj = [zeros(E*N*F+D*D,1); zeros(P,1); zeros(tau_len,1); zeros(kappa_len, 1); zeros(theta_len, 1)];

% dimension checks
num_rows = offset + 2*(tau_len+kappa_len+theta_len);
num_cols = E*N*F + D*D + P + tau_len + kappa_len + theta_len;
assert(size(A,1)==num_rows);
assert(size(A,2)==num_cols); assert(length(LB)==num_cols); assert(length(UB)==num_cols); assert(length(obj)==num_cols);


%%%%% quadratic term of the objective function %%%%%
disp('generating Q');
idx_loss = E*N*F+D*D+1:E*N*F+D*D+P;
idx_tau = E*N*F+D*D+P+1:E*N*F+D*D+P+tau_len;
idx_kappa = E*N*F+D*D+P+tau_len+1:E*N*F+D*D+P+tau_len+kappa_len;
idx_theta = E*N*F+D*D+P+tau_len+kappa_len+1:E*N*F+P+D*D+tau_len+kappa_len+theta_len;

Qi = [idx_loss, idx_tau, idx_kappa, idx_theta];
Qj = [idx_loss, idx_tau, idx_kappa, idx_theta];
Qv = [ones(1,P), CFG.C2.tau * ones(1,tau_len), CFG.C2.kappa * ones(1,kappa_len), CFG.C2.theta * ones(1,theta_len)];
Q = sparse(Qi, Qj, Qv, length(obj), length(obj));
clear Qi Qj Qv;

sA = whos('A');
sQ = whos('Q');
total_size = sA.bytes+sQ.bytes+numel(LB)*3*8 + numel(b)*8;
fprintf('total problem size: %1.2f mb\n', total_size/1024^2);

used_mem = whos;
fprintf('used memory: %1.2f mb\n', sum([used_mem.bytes])/1024^2);


%%%%% solve optimisation problem %%%%%
optimizer = 'cplex';
switch optimizer
%switch CFG.optimizer
 case 'cplex'
  lpenv = cplex_license(1,1);
  %how = lp_set_param(lpenv,'CPX_PARAM_PREDUAL', 0, 1); 
  how = lp_set_param(lpenv, 'CPX_PARAM_AGGFILL', 10000, 1); 
  how = lp_set_param(lpenv, 'CPX_PARAM_AGGIND', 10000, 1); 
  how = lp_set_param(lpenv, 'CPX_PARAM_PREPASS', 10000, 1); 
  how = lp_set_param(lpenv, 'CPX_PARAM_BARCOLNZ', 10, 1); 
  how = lp_set_param(lpenv, 'CPX_PARAM_BAREPCOMP', 1e-6, 1); 
  how = lp_set_param(lpenv, 'CPX_PARAM_PRECOMPRESS', 1, 1);  
  %how = lp_set_param(lpenv, 'CPX_PARAM_BARTHREADS', 1, 1); 
  %how = lp_set_param(lpenv, 'CPX_PARAM_THREADS', 4, 1); 
  %how = lp_set_param(lpenv, 'CPX_PARAM_THREADS', 8, 1); 
  if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
  [xopt, lambda, how] = qp_solve(lpenv, Q, obj, sparse(A), b, LB, UB, num_rows - 2*(tau_len + kappa_len + theta_len), double(CFG.VERBOSE>0), 'bar');
  if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end
  if ~isequal(how,'OK')
    warning('CPLEX: %s\n', how);
  end
  [lpenv, status] = cplex_quit(lpenv,0);
  lpenv = 0;
 case 'mosek'
  idx_eq = 1:num_rows-2*(tau_len + kappa_len + theta_len);
  idx_neq = num_rows-2*(tau_len + kappa_len + theta_len)+1:num_rows;
  if CFG.VERBOSE>0,
    display_mode = 'iter';
  else
    display_mode = 'off';
  end
  if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
  [xopt, lambda] = quadprog(Q, obj, sparse(A(idx_neq,:)), b(idx_neq), sparse(A(idx_eq,:)), b(idx_eq), LB, UB, [], optimset('Display', display_mode));
  if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end
 otherwise
  error('unknown optimizer %s', CFG.optimizer);
end

weights = reshape(xopt(1:E*N*F)', F, N, E);
dists = reshape(xopt(E*N*F+1:E*N*F+D*D), D, D);
xis = xopt(E*N*F+D*D+1:E*N*F+D*D+P);
taus = xopt(idx_tau);
kappas = xopt(idx_kappa);
thetas = xopt(idx_theta);

objective.tau = 0.5*CFG.C2.tau*sum(taus.^2);
objective.kappa = 0.5*CFG.C2.kappa*sum(kappas.^2);
objective.theta = 0.5*CFG.C2.theta*sum(thetas.^2);
objective.all = sum(xopt.*obj) + 0.5*xopt'*Q*xopt;

objective

