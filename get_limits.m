function limits = get_limits(max_lmt, num_lmt)
% limits = get_limits(max_lmt, num_lmt)
%
% -- input --
% max_lmt: maximal supporting point of PLiF
% num_lmt: number of supporting points
%
% -- ouput --
% limits: supporting points of PLiF


limits = round(linspace(0, sqrt(max_lmt), num_lmt).^2);
limits(1) = 1;
limits(end) = inf;