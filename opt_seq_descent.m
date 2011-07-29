function [seq_weights, obj, fval] = opt_seq_descent(CFG, seq_weights, seq_feat, exon_mask, tmp_seq_weights, weights, coverage, R_const)
% OPT_SEQ_DESCENT   Determines the optimal sequence weights.
%
%   [seq_weights, obj, fval] = opt_seq_descent(CFG, seq_weights, seq_feat, exon_mask, tmp_seq_weights, weights, coverage, R_const)
%
%   -- input --
%   CFG:             configuration struct
%   seq_weights:     weights for sequence normalisation
%   seq_feat:        vector of encoded sequence features
%   exon_mask:       P x T matrix of profile correction
%   tmp_seq_weights: S*T x T matrix of sequence weights for each transcript
%   weights:         weights of transcripts
%   coverage:        vector of observed exon coverage
%   R_const:         constant residue
%
%   -- output --
%   seq_weights:     weights for sequence normalisation
%   obj:             objective value at optimum
%   fval:            objective value at each step
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


S = 0;
for o = 1:CFG.seq.order,
  S = S + (CFG.seq.half_win_size*2-o+1) * 4^o;
end
T = length(weights);
ts_idx = zeros(1,S*T);
for n = 1:S,
  ts_idx([1:S:length(ts_idx)]+n-1) = n+[0:(S*T+S):S*T*T];
end

max_iter = 1;

fval = 1e100*ones(1, S);
fval_old = zeros(1, S);
seq_weights_old = zeros(size(seq_weights));

if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
iter = 1;
while 1
  cnt = 1;
  for k = 1:S,
     %%% Rbe: residue for beta
     Rbe = exon_mask.*seq_feat(:,k:S:S*T)*weights';
     %%% R1: residue for beta independent variables
     idx_wo_be = setdiff([1:S*T], [k:S:S*T]);
     tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
     R1 = exon_mask.*(seq_feat(:,idx_wo_be)*tmp_seq_weights(idx_wo_be,:))*weights' - ...
          coverage;
     %%% S1: residue of quadratic term
     S1 = sum(Rbe.^2);
     %%% S2: residue of linear term
     S2 = 2*sum(Rbe'*R1);
     %%% S3: constant term
     S3 = sum(R1.^2) + R_const;
     %%% calculation and clipping of beta
     be_new = -0.5*S2/S1;
     if be_new < 0
       seq_weights(k) = 0.0;
     else
       seq_weights(k) = be_new;
     end
     fval(cnt) = quad_fun(seq_weights(k), S1, S2, S3);
     tmp_seq_weights(ts_idx) = repmat(seq_weights, T, 1);
     obj_alt = sum((exon_mask.*(seq_feat*tmp_seq_weights)*weights'-coverage).^2) + R_const;
     assert(abs(fval(cnt)-obj_alt)>-1e-3); % objective should be indentical to not-expanded objective
     cnt = cnt + 1;
  end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\n', iter, fval(end), norm(seq_weights_old-seq_weights)); end
  if norm(fval_old-fval)<1e-5 || norm(seq_weights_old-seq_weights)<1e-5 || iter>=max_iter,
    break;
  end
  iter = iter + 1;
end
assert(all(fval(1:end-1)-fval(2:end)>-1e-3)); % objective should decrease at every step
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

obj = fval(end);

