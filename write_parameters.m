function write_parameters(PAR, fname)
% WRITE_PARAMETERS   Writes parameters to file
%
%   write_parameters(PAR, fname)
%
%   -- input --
%   PAR:   parameter struct
%   fname: name of file
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


[fd msg] = fopen(fname, 'w+');
assert(fd~=-1);
fprintf(fd, 'C_I\t%i\n', PAR.CFG.C_I);
fprintf(fd, 'C_F\t%i\n', PAR.CFG.C_F);
fprintf(fd, 'C_N\t%i\n', PAR.CFG.C_N);
fprintf(fd, 'learn_profiles\t%i\n', PAR.learn_profiles);
fprintf(fd, 'load_profiles\t%i\n', PAR.load_profiles);
fprintf(fd, 'profiles\t%s\n', PAR.profiles_fn);
fclose(fd);