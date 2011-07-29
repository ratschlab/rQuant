function plot_sequence_weights(weights, order, win_size)
% PLOT_SEQUENCE_WEIGHTS   Visualises sequence weights.
%
%   plot_sequence_weights(weights, order, win_size)
%
%   -- input --
%   weights:  estimated sequence weights
%   order:    order of k-mers
%   win_size: size of sequence window
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


acgt = 'ACGT';
num_nt = length(acgt);

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

offset = 0;
for o = 1:order,
  offset(o+1) = offset(o) + (win_size-o+1) * 4^o;
  figure(o);
  imagesc(reshape(weights(offset(o)+1:offset(o+1)), 4^o, win_size-o+1));
  xlabel('position'); ylabel('k-mer')
  ylim([0.5 4^o+0.5]); set(gca, 'YTick', 1:1:4^o);
  set(gca, 'YTickLabel', nucleotides{o});
  colorbar;
end