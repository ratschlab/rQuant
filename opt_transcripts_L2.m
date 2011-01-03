function [weights, obj, fval] = opt_transcripts_L2(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, max_iter, weights0)
% [weights, obj, fval] = opt_transcripts_L2(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, max_iter)
%
% -- input --
% CFG: configuration struct
% coverage: vector of observed exon coverage 
% exon_mask: PxT matrix modelling 'importance' of a position within a transcript
% intron_count: vector of observed intron confirmation
% intron_mask: IxT matrix defining whether an intron belongs to a particular transcript
% C_w: regularisation parameter per transcript (T x 1)
% max_iter: maximal number of iterations
%
% -- output --
% weights: weights of transcripts
% obj: objective value at optimum
% fval: 

if nargin<7
  max_iter = 1e100;
end

T = size(exon_mask,2); % number of transcripts
I = size(intron_mask,1); % number of introns

exon_count = sum(coverage,2);

if nargin<8
  weights = full(mean(coverage)/T*ones(1,T));
else
  weights = weights0;
end
%weights = rand(1,T);
weights_old = zeros(1,T);
fval = 1e100*ones(1,T);
fval_old = zeros(1,T);

LB = 0;
UB = full(mean(coverage));

if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
iter = 1;
if T==1
  RE = -exon_count;
  RI = -intron_count;
  S1 = sum(exon_mask.^2) + C_w; if I>0, S1 = S1 + sum(intron_mask.^2); end
  S2 = 2*sum(exon_mask'*RE); if I>0, S2 = S2+2*sum(intron_mask'*RI); end
  S3 = sum(RE.^2); if I>0, S3 = S3 + sum(RI.^2); end
  weights = -0.5*S2/S1;
  fval = quad_fun(weights(t), S1, S2, S3);
  %[weights fval exitflag] = fminbnd(@(w) quad_fun(w, S1, S2, S3), LB, UB);
else
  while 1
    weights_old = weights;
    fval_old = fval;
    for t = 1:T,
      idx_wo_t = setdiff(1:T,t);
      RE = exon_mask(:,idx_wo_t)*weights(idx_wo_t)'-exon_count;
      if I>0, RI = intron_mask(:,idx_wo_t)*weights(idx_wo_t)'-intron_count; end
      S1 = sum(exon_mask(:,t).^2) + C_w(t); if I>0, S1 = S1 + sum(intron_mask(:,t).^2); end
      S2 = 2*sum(exon_mask(:,t)'*RE); if I>0, S2 = S2+2*sum(intron_mask(:,t)'*RI); end
      S3 = sum(RE.^2) + weights(idx_wo_t).^2*C_w(idx_wo_t); if I>0, S3 = S3 + sum(RI.^2); end
      weights(t) = -0.5*S2/S1;
      fval(t) = quad_fun(weights(t), S1, S2, S3);
      %[weights(t) fval(t) exitflag] = fminbnd(@(w) quad_fun(w, S1, S2, S3), LB, UB);
    end
    if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\n', iter, fval(end), norm(weights_old-weights)); end
    if norm(fval_old-fval)<1e-5 | norm(weights_old-weights)<1e-5 | iter>=max_iter,
      break;
    end
    iter = iter + 1;
  end
end
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

obj = fval(end);