function [weights, obj, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, R_const, max_iter, weights0, reg)
% OPT_TRANSCRIPTS_DESCENT   Determines the optimal transcript weights.
%
%   [weights, obj, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, R_const, max_iter, weights0, reg)
%
%   -- input --
%   CFG:          configuration struct
%   coverage:     vector of observed exon coverage 
%   exon_mask:    PxT matrix modelling 'importance' of a position within a transcript
%   intron_count: vector of observed intron confirmation
%   intron_mask:  IxT matrix defining whether an intron belongs to a particular transcript
%   C_w:          regularisation parameter per transcript (T x 1)
%   R_const:      constant residue
%   max_iter:     maximal number of iterations (optional)
%   weights0:     initialisation values of the weights (optional)
%   reg:          regularisation method (optional)
%
%   -- output --
%   weights:     weights of transcripts
%   obj:         objective value at optimum
%   fval:        objective value at each step
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%


T = size(exon_mask,2); % number of transcripts
I = size(intron_mask,1); % number of introns

exon_count = sum(coverage,2);

if nargin<7
  R_const = 0;
end

if nargin<8
  max_iter = 1e100;
end

if nargin<9
  weights = full(mean(coverage)/T*ones(1,T));
else
  weights = weights0;
end
weights_old = zeros(1,T);
fval = 1e100*ones(1,T);
fval_old = zeros(1,T);

if nargin<10
  reg = 'L1';
end

LB = 0.0;
UB = full(mean(coverage));

cnt = 0;
if CFG.VERBOSE>0, fprintf(1, '\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
if T==1
  RE = -exon_count;
  RI = -intron_count;
  switch reg
   case 'L1'
    S1 = sum(exon_mask(:,1).^2); if I>0, S1 = S1 + CFG.C_I*sum(intron_mask(:,1).^2); end
    S2 = 2*sum(exon_mask(:,1)'*RE) + C_w(1); if I>0, S2 = S2+2*CFG.C_I*sum(intron_mask(:,1)'*RI); end
    S3 = sum(RE.^2) + R_const; if I>0, S3 = S3 + CFG.C_I*sum(RI.^2); end
   case 'L2'
    S1 = sum(exon_mask(:,1).^2) + C_w(1); if I>0, S1 = S1 + CFG.C_I*sum(intron_mask(:,1).^2); end
    S2 = 2*sum(exon_mask(:,1)'*RE); if I>0, S2 = S2+2*CFG.C_I*sum(intron_mask(:,1)'*RI); end
    S3 = sum(RE.^2) + R_const; if I>0, S3 = S3 + CFG.C_I*sum(RI.^2); end
  end
  w_new = -0.5*S2/S1;
  % clipping of w_t
  if w_new < 0
    weights(1) = LB;
  else
    weights(1) = w_new;
  end
  fval(1) = quad_fun(weights(1), S1, S2, S3);
  obj_alt = sum((exon_mask*weights'-coverage).^2) + R_const;
  if I>0, obj_alt = obj_alt + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
  switch reg
   case 'L1'
    obj_alt = obj_alt + abs(weights*C_w);
   case 'L2'
    obj_alt = obj_alt + weights.^2*C_w;
  end
  if ~(abs(fval(1)-obj_alt)<1e-3) % objective should be indentical to not-expanded objective
    cnt = cnt + 1;
    if CFG.VERBOSE>1, fprintf(1, 'objectives differ %.6f (tscp %i)\n', abs(fval(t)-obj_alt), t); end
  end
else
  iter = 1;
  while 1
    weights_old = weights;
    fval_old = fval;
    for t = 1:T,
      idx_wo_t = setdiff(1:T,t);
      RE = exon_mask(:,idx_wo_t)*weights(idx_wo_t)'-exon_count;
      if I>0, RI = intron_mask(:,idx_wo_t)*weights(idx_wo_t)'-intron_count; end
      switch reg
       case 'L1'
        S1 = sum(exon_mask(:,t).^2); if I>0, S1 = S1 + CFG.C_I*sum(intron_mask(:,t).^2); end
        S2 = 2*sum(exon_mask(:,t)'*RE) + C_w(t); if I>0, S2 = S2+2*CFG.C_I*sum(intron_mask(:,t)'*RI); end
        S3 = sum(RE.^2) + R_const + abs(weights(idx_wo_t))*C_w(idx_wo_t); if I>0, S3 = S3 + CFG.C_I*sum(RI.^2); end
       case 'L2'
        S1 = sum(exon_mask(:,t).^2) + C_w(t); if I>0, S1 = S1 + CFG.C_I*sum(intron_mask(:,t).^2); end
        S2 = 2*sum(exon_mask(:,t)'*RE); if I>0, S2 = S2+2*CFG.C_I*sum(intron_mask(:,t)'*RI); end
        S3 = sum(RE.^2) + R_const + weights(idx_wo_t).^2*C_w(idx_wo_t); if I>0, S3 = S3 + CFG.C_I*sum(RI.^2); end
      end
      w_new = -0.5*S2/S1;
      % clipping of w_t
      if w_new < 0
        weights(t) = LB;
      else
        weights(t) = w_new;
      end
      fval(t) = quad_fun(weights(t), S1, S2, S3);
      obj_alt = sum((exon_mask*weights'-coverage).^2) + R_const;
      if I>0, obj_alt = obj_alt + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
      switch reg
       case 'L1'
        obj_alt = obj_alt + abs(weights*C_w);
       case 'L2'
        obj_alt = obj_alt + weights.^2*C_w;
      end
      if ~(abs(fval(t)-obj_alt)<1e-3) % objective should be indentical to not-expanded objective
        cnt = cnt + 1;
        if CFG.VERBOSE>1, fprintf(1, 'objectives differ %.6f (tscp %i)\n', abs(fval(t)-obj_alt), t); end
      end
    end
    if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\n', iter, fval(end), norm(weights_old-weights)); end
    if norm(fval_old-fval)<1e-5 || norm(weights_old-weights)<1e-5 || iter>=max_iter,
      break;
    end
    iter = iter + 1;
  end
end
if CFG.VERBOSE>0 && cnt>0, fprintf(1, 'objectives differ for %i transcripts\n', cnt); end
if 0%mean(fval(1:end-1)-fval(2:end)>-1e-2)<0.95
  fprintf(1, 'objective non-decreasing\n');
end
assert(all(fval(1:end-1)-fval(2:end)>-1e-2));
if CFG.VERBOSE>0, fprintf(1, 'Took %.1fs.\n', toc); end

obj = fval(end);