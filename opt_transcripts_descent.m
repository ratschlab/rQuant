function [weights, obj, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, max_iter, weights0, reg)
% [weights, obj, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, C_w, max_iter)
%
% -- input --
% CFG: configuration struct
% coverage: vector of observed exon coverage 
% exon_mask: PxT matrix modelling 'importance' of a position within a transcript
% intron_count: vector of observed intron confirmation
% intron_mask: IxT matrix defining whether an intron belongs to a particular transcript
% C_w: regularisation parameter per transcript (T x 1)
% max_iter: maximal number of iterations (optional)
% weights0: initialisation values of the weights (optional)
% reg: regularisation method (optional)
%
% -- output --
% weights: weights of transcripts
% obj: objective value at optimum
% fval: objective value at each step

if nargin<9
  reg = 'L1';
end
  
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
weights_old = zeros(1,T);
fval = 1e100*ones(1,T);
fval_old = zeros(1,T);

LB = 0.0;
UB = full(mean(coverage));

if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
if T==1
  RE = -exon_count;
  RI = -intron_count;
  switch reg
   case 'L1'
    S1 = sum(exon_mask(:,1).^2); if I>0, S1 = S1 + sum(intron_mask(:,1).^2); end
    S2 = 2*sum(exon_mask(:,1)'*RE) + C_w(1); if I>0, S2 = S2+2*sum(intron_mask(:,1)'*RI); end
    S3 = sum(RE.^2); if I>0, S3 = S3 + sum(RI.^2); end
   case 'L2'
    S1 = sum(exon_mask(:,1).^2) + C_w(1); if I>0, S1 = S1 + sum(intron_mask(:,1).^2); end
    S2 = 2*sum(exon_mask(:,1)'*RE); if I>0, S2 = S2+2*sum(intron_mask(:,1)'*RI); end
    S3 = sum(RE.^2); if I>0, S3 = S3 + sum(RI.^2); end
  end
  w_new = -0.5*S2/S1;
  % clipping of w_t
  if w_new < 0
    weights(1) = LB;
  else
    weights(1) = w_new;
  end
  fval(1) = quad_fun(weights(1), S1, S2, S3);
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
        S1 = sum(exon_mask(:,t).^2); if I>0, S1 = S1 + sum(intron_mask(:,t).^2); end
        S2 = 2*sum(exon_mask(:,t)'*RE) + C_w(t); if I>0, S2 = S2+2*sum(intron_mask(:,t)'*RI); end
        S3 = sum(RE.^2) + abs(weights(idx_wo_t))*C_w(idx_wo_t); if I>0, S3 = S3 + sum(RI.^2); end
       case 'L2'
        S1 = sum(exon_mask(:,t).^2) + C_w(t); if I>0, S1 = S1 + sum(intron_mask(:,t).^2); end
        S2 = 2*sum(exon_mask(:,t)'*RE); if I>0, S2 = S2+2*sum(intron_mask(:,t)'*RI); end
        S3 = sum(RE.^2) + weights(idx_wo_t).^2*C_w(idx_wo_t); if I>0, S3 = S3 + sum(RI.^2); end
      end
      w_new = -0.5*S2/S1;
      % clipping of w_t
      if w_new < 0
        weights(t) = LB;
      else
        weights(t) = w_new;
      end
      fval(t) = quad_fun(weights(t), S1, S2, S3);
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