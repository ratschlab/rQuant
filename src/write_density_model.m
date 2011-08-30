function write_density_model(profile_weights profiles_fn_out)
% WRITE_DENSITY_MODEL   Writes density models from matrix into text file.
%
%   write_density_model(profile_weights profiles_fn_out)
%
%   -- input --
%   profile_weights: learned weights of profile functions
%   profiles_fn_out: name of output file that stores profiles
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
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