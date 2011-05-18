function write_density_model(profile_weights, intron_dists, profiles_fn_out, intron_dists_fn)
% WRITE_DENSITY_MODEL   Writes density models from matrices into text files.
%
%   write_density_model(profile_weights, intron_dists, profiles_fn_out, intron_dists_fn)
%
%   -- input --
%   profile_weights: learned weights of profile functions
%   intron_dists:    learned distances to closest intron
%   profiles_fn_out: name of output file that stores profiles
%   intron_dists_fn: name of output files that stores intron distances
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


% write profile_weights to output file
[fd msg] = fopen(profiles_fn_out, 'w+');
assert(fd~=-1);
for n = 1:size(profile_weights,1)
  for m = 1:size(profile_weights,2),
    fprintf(fd, '%4.8f\t', profile_weights(n,m));
  end
  fprintf(fd, '\n');
end
fclose(fd);

% write intron dists to output file
[fd msg] = fopen(intron_dists_fn, 'w+');
assert(fd~=-1);
for n = 1:size(intron_dists,1)
  for m = 1:size(intron_dists,2),
    fprintf(fd, '%4.8f\t', intron_dists(n,m));
  end
  fprintf(fd, '\n');
end
fclose(fd);

% plot profile weights
%figure(); plot(profile_weights);
%print('-dpng', '-r150', strrep(profiles_fn_out, '.txt', 'png'));