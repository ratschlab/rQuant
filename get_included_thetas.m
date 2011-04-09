function pw_nnz = get_included_thetas(CFG)
% pw_nnz = get_included_thetas(CFG)
%
% -- input --
% CFG: configuration struct
%
% -- ouput --
% pw_nnz: F x N logical matrix indicating included supporting points (thetas)

F = CFG.num_plifs;
N = size(CFG.transcript_len_ranges,1);
% supporting points of PLiFs
lmt = get_limits(CFG.max_side_len, F/2);
pw_nnz = false(F, N);
for n = 1:N,
  fidx = find(CFG.transcript_len_ranges(n,2)/2>lmt, 1, 'last');
  pw_nnz([1:fidx+1, F-fidx:F],n) = true;
end