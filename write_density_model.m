%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function write_density_model(profile_weights, intron_dists, profiles_fn_out, intron_dists_fn)
% write_density_model(profile_weights, intron_dists, profiles_fn_out, intron_dists_fn)
%
% -- input --
% profile_weights: learned weights of profile functions
% intron_dists: learned distances to closest intron
% profiles_fn_out
% intron_dists_fn
  
% write profile_weights to output file
[fd msg] = fopen(profiles_fn_out, 'w+');
disp(msg);
for n = 1:size(profile_weights,1)
  for m = 1:size(profile_weights,2),
    fprintf(fd, '%4.8f\t', profile_weights(n,m));
  end
  fprintf(fd, '\n');
end
fclose(fd);

% write intron dists to output file
[fd msg] = fopen(intron_dists_fn, 'w+');
disp(msg);
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