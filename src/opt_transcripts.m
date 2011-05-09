function [weights, betas, xis, loss] = opt_transcripts(CFG, gene, coverage, exon_mask, excluded_reads, intron_count, intron_mask, lpenv)
% OPT_TRANSCRIPTS   Determines the optimal transcript weights.
%
%   [weights, betas, xis, loss] = opt_transcripts(CFG, gene, coverage, exon_mask, excluded_reads, intron_count, intron_mask, lpenv)
%
%   -- input --
%   CFG:          configuration struct
%   gene:         struct defining a gene with start, stops, exons etc.
%   coverage:     vector of observed exon coverage 
%   exon_mask:    PxT matrix modelling 'importance' of a position within a transcript
%   intron_count: vector of observed intron confirmation
%   intron_mask:  IxT matrix defining whether an intron belongs to a particular transcript
%   lpenv:        license for optimiser
%
%   -- output --
%   weights:      weights of transcripts
%   betas:        weights of reads
%   xis:          slack variables
%   loss:         exon, intron and total loss
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


%%%%% define optimisation problem %%%%%
% min_{x=[w,b,xis]}  sum_{t=1}^{T} C_1*w_{t}
%                  + sum_{t=1}^{T} sum_{r=1}^{R} C_beta b_{t,r}
%                  + sum_{p=1}^{P} sum_{t=1}^{T} xi_{t,p}
%                  + sum_{i=1}^{I} sum_{i in I} xi_{i}
% s.t.
% for all p in P,
%     all t in T: -em_{t,p}*w_{t} + sum_{r=1}^{R} rm_{p,t,r}*b_{t,r} - xi_{t,p} = 0
% for all r in R:                              sum_{t=1}^{T} b_{t,r}            = 1
% for all i in I: -sum_{t=1}^{T} im_{t,i}*w_{t}                      - xi_{i}   = -ic_{i}
% 
% with bounds:
%    0 <= w_{t}    <= INF     for all t=1..T
%    0 <= b_{t,r}  <= 1       for all t=1..T, r=1..R
% -INF <= xi_{p}   <= INF     for all p=1..P
% -INF <= xi_{i}   <= INF     for all i=1..I
%
% T: number of transcripts
% P: number of positions
% I: number of introns
% R: number of reads


INF = 1e20;

% subsample if problem gets to large
max_num_exm = 20000;
if sum(max(exon_mask,[],2))>max_num_exm
  p = max_num_exm/sum(max(exon_mask,[],2));
  if CFG.VERBOSE>0, fprintf('subsampling positions p=%1.2f\n', p); end
  pidx = randperm(size(exon_mask,1));
  pidx = sort(pidx(1:round(p*size(exon_mask,1))));
  exon_mask = exon_mask(pidx,:);
  coverage = coverage(pidx, :);
  assert(size(exon_mask,1)>0 & size(exon_mask,2)>0);
end
if size(coverage,2)>max_num_exm/2
  if CFG.VERBOSE>0, fprintf('summing up %i reads\n', size(coverage,2)); end
  coverage = sum(coverage,2);
  for t = 1:length(gene.transcripts),
    excluded_reads{t} = [];
  end
end


T = size(exon_mask,2);   % number of transcripts
P = size(exon_mask,1);   % number of positions
I = size(intron_mask,1); % number of introns
R = size(coverage,2);    % number of reads

size_sparse = P*(T+T*R+T) + T*R;

b_exon = zeros(P*T+R,1);
Ac = 0; % matrix counter 

if CFG.VERBOSE>1, fprintf(1, 'number of nz in A: %i\n', size_sparse); end
Ai = zeros(1,size_sparse); % row indices
Aj = zeros(1,size_sparse); % column indices
Av = zeros(1,size_sparse); % matrix values

% exon part
% T*P constraints
for p = 1:P,
  Ai(Ac+1:Ac+T)   = T*(p-1)+1:T*(p-1)+T;
  Aj(Ac+1:Ac+T)   = 1:T;
  Av(Ac+1:Ac+T)   = -exon_mask(p,:); 
  Ac              = Ac + T;
  Ai(Ac+1:Ac+T*R) = reshape(repmat([T*(p-1)+1:T*(p-1)+T]',1,R)', R*T, 1)';
  Aj(Ac+1:Ac+T*R) = T + [1:T*R];
  Av(Ac+1:Ac+T*R) = repmat(coverage(p,:), 1, T);
  Ac              = Ac + T*R;
  Ai(Ac+[1:T])    = T*(p-1)+1:T*(p-1)+T;
  Aj(Ac+[1:T])    = T + T*R + [T*(p-1)+1:T*(p-1)+T];
  Av(Ac+[1:T])    = -1;
  Ac              = Ac + T;
  b_exon(T*(p-1)+1:T*(p-1)+T,1) = 0;
end
% sum_{t} b_{t,r} = 1 
for r = 1:R,
  Ai(Ac+[1:T]) = P*T+r;
  Aj(Ac+[1:T]) = [T+r:R:T+T*R];
  Av(Ac+[1:T]) = 1; 
  b_exon(P*T+r,1) = 1;
  Ac = Ac+T;
end
assert(Ac==size_sparse);
idx = find(Av);
A_exon = sparse(Ai(idx), Aj(idx), Av(idx), P*T+R, T+T*R+P*T+I);
% intron part 
b_intron = zeros(I, 1);
Ac = 0; % matrix counter 
Ai = zeros(1,I*(T+1)); % row indices
Aj = zeros(1,I*(T+1)); % column indices
Av = zeros(1,I*(T+1)); % matrix values
for idx = 1:I,
  Ai(Ac+1:Ac+T)   = idx;
  Aj(Ac+1:Ac+T)   = 1:T;
  Av(Ac+1:Ac+T)   = -intron_mask(idx,:); 
  Ac              = Ac + T;
  Ai(Ac+1)        = idx;
  Aj(Ac+1)        = T + T*R + P*T + idx;
  Av(Ac+1)        = -1;
  Ac              = Ac + 1;
  b_intron(idx,1) = -intron_count(idx);
end
A_intron = sparse(Ai, Aj, Av, I, T+T*R+P*T+I);

A = [A_exon; A_intron]; % (P*T+R+I)x(T+T*R+P*T+I) matrix
b = [b_exon; b_intron]; % (P*T+R+I) vector


%%%%% bounds %%%%%
LB = [zeros(T+T*R,1); -INF*ones(P*T+I,1)]; % lower bounds for xopt
UB = [INF*ones(T,1); ones(T*R,1); INF*ones(P*T+I,1)]; % upper bounds for xopt

C_betas = zeros(R,T);
for t = 1:T,
  C_betas(excluded_reads{t},t) = 10^6;
end
C_betas = reshape(C_betas, T*R, 1);


%%%%% linear term of the objective function %%%%%
if ~isfield(gene, 'C1'), gene.C1(1:length(gene.transcripts)) = 1; end
if ~isfield(gene, 'transcript_len'), gene.transcript_len = gene.transcript_length; end
obj = [(CFG.C1*gene.transcript_len.*(gene.C1.^2))'; C_betas; zeros(P*T+I,1)];

% dimension checks
num_rows = P*T+R+I;
num_cols = T+T*R+P*T+I;
assert(size(A,1)==num_rows);
assert(size(A,2)==num_cols); assert(length(LB)==num_cols); assert(length(UB)==num_cols); assert(length(obj)==num_cols);


%%%%% quadratic term of the objective function %%%%%
intron_factor = (10*CFG.read_len/2);
Q = spzeros(length(obj));
for idx = T+T*R+1:T+T*R+T*P,
  Q(idx,idx) = 1;
end
for idx = T+T*R+T*P+1:T+T*R+T*P+I,
  Q(idx,idx) = intron_factor;
end


%%%%% solve optimisation problem %%%%%
switch CFG.optimizer
 case 'mosek'
  idx_neq = [];
  idx_eq = 1:num_rows;
  if CFG.VERBOSE>1,
    display_mode = 'iter';
  else
    display_mode = 'off';
  end
  if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
    switch determine_engine
     case 'octave'
      xopt = mosek_qp([], Q, obj, sparse(A(idx_eq,:)), b(idx_eq), LB, UB, -inf*ones(size(b(idx_neq))), sparse(A(idx_neq,:)), b(idx_neq), 4, CFG.VERBOSE-1);
     case 'matlab'
      [xopt, lambda] = quadprog(Q, obj, sparse(A(idx_neq,:)), b(idx_neq), sparse(A(idx_eq,:)), b(idx_eq), LB, UB, [], optimset('Display', display_mode));
     otherwise
      error('unknown engine %s', determine_engine);  
    end
    if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end
 otherwise
  error('unknown optimizer %s', CFG.optimizer);
end
  
loss.exons = 0.5*sum(xopt(T+T*R+1:T+T*R+P*T).^2);
loss.introns = 0.5*intron_factor*sum(xopt(T+T*R+P*T+1:T+T*R+P*T+I).^2);
loss.all = sum(xopt.*obj) + loss.exons + loss.introns;

weights = xopt(1:T)';
betas = reshape(xopt(T+1:T+T*R), R, T)';
xis = xopt(T+T*R+1:end);