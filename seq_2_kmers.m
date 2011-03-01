%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2007-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2007-2010 Max Planck Society
%

function kmers = seq_2_kmers(seq, order, win_size)
% kmers = seq_2_kmers(seq, order, win_size)
%
% -- input --
% seq: sequence
% order: maximal substring length
% win_size: size of subsequence window
%
% -- output --
% kmers: matrix of positional substring occurrence


assert(all(all(seq=='A' | seq=='C' | seq=='G' | seq=='T')));
s_len = length(seq);
acgt = 'ACGT';
num_nt = length(acgt);

% determine number of k-mer occurences in whole sequence
offset_all = zeros(1, order+1);
for o = 1:order,
  offset_all(o+1) = offset_all(o) + num_nt^o;
end
kmer_dim = offset_all(end);

kmers_all = zeros(kmer_dim, s_len);
for n = 1:num_nt,
  kmers_all(n,:) = double(seq==acgt(n));
end  
for o = 2:order,
  for n = offset_all(o-1)+1:offset_all(o),
    kmers_all(offset_all(o)+(n-offset_all(o-1)-1)*num_nt+[1:num_nt],1:end-o+1) = repmat(kmers_all(n,1:end-o+1), num_nt, 1) .* kmers_all(1:num_nt,o:end);
  end
end

% save number of k-mer occurences per subsequence for each position
offset = zeros(1, order+1);
for o = 1:order,
  offset(o+1) = offset(o) + (win_size-o+1) * num_nt^o;
end
kmer_dim = offset(end);
kmers = zeros(kmer_dim, length(seq)-win_size);
for n = win_size/2+1:length(seq)-win_size/2,
  for o = 1:order,
    kmers(offset(o)+1:offset(o+1),n-win_size/2) = reshape(kmers_all(offset_all(o)+1:offset_all(o+1), n-win_size/2:n+win_size/2-1-o+1), (win_size-o+1) * num_nt^o, 1);
  end
end

exp_sum = 0;
for o = 1:order,
  exp_sum = exp_sum + win_size-o+1;
end
assert(all(sum(kmers) == exp_sum));
