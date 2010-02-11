function kmers = seq_2_kmers(seq, order)
% kmers = seq_2_kmers(seq, order)
%
% -- input --
% seq: matrix of sequences
% order: maximal substring length
%
% -- output --
% kmers: matrix of positional substring occurrence


assert(all(all(seq=='A' | seq=='C' | seq=='G' | seq=='T')));
s_len = size(seq,1);
offset = 0;
acgt = 'ACGT';
num_nt = length(acgt);
for o = 1:order,
  offset(o+1) = offset(o) + (s_len-o+1) * num_nt^o;
end
kmer_dim = offset(end);

nucleotides{1} = acgt';
for o = 2:order,
  nucleotides{o} = char(zeros(length(nucleotides{o-1}) * num_nt, o));
  for s = 1:length(nucleotides{o-1}),
    for n = 1:num_nt,
      row = (s-1) * num_nt + n;
      nucleotides{o}(row,:) = [nucleotides{o-1}(s,:) acgt(n)];
    end
  end
end

kmers = zeros(kmer_dim, size(seq,2));
for o = 1:order,
  for s = 1:s_len-o+1,
    for n = 1:length(nucleotides{o});
      idx = strmatch(nucleotides{o}(n,:), seq(s:s+o-1,:)');
      row = offset(o) + (s-1)*num_nt^o + n;
      kmers(row,idx) = 1;
    end
  end
end

exp_sum = 0;
for o = 1:order,
  exp_sum = exp_sum + s_len-o+1;
end
assert(all(sum(kmers) == exp_sum));


