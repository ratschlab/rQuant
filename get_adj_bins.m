function p_adj = get_adj_bins(CFG)
% p_adj = get_adj_bins(CFG)
% 
% -- input --
% CFG: configuration struct
%
% -- output --
% p_adj: F*N x 2 matrix of adjacent supporting points (1: prevoius, 2: next)


pw_nnz = get_included_thetas(CFG);
F = CFG.num_plifs;
N = size(CFG.transcript_len_ranges,1);

p_adj = zeros(2, F*N);
p_adj(1,:) = 0:F*N-1;
p_adj(1,1+F*[0:N-1]) = p_adj(1,2+F*[0:N-1]);
fidx = find((~pw_nnz))+1;
for f = 1:length(fidx),
  p_adj(1,fidx(f)) = p_adj(1,fidx(f)-1);
end
p_adj(2,:) = 2:F*N+1;
p_adj(2,F+F*[0:N-1]) = p_adj(2,F-1+F*[0:N-1]);
fidx = find(~pw_nnz)-1;
for f = length(fidx):-1:1,
  p_adj(2,fidx(f)) = p_adj(2,fidx(f)+1);
end