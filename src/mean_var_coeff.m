function [Q1 Q2] = mean_var_coeff(pred, eval)
% MEAN_VAR_COEFF   Calculates the ratio of mean to variance.
%
%   [Q1 Q2] = mean_var_coeff(pred, eval)
%
%   -- input --
%   pred: predicted values
%   eval: actual values
%
%   -- output --
%   Q1:   absolute variability coefficient
%   Q2:   squared variability coefficient
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


Q1 = mean(abs(pred - eval)) / mean(abs(eval - mean(eval)));
%Q1 = mean(abs(pred - eval)) / mean(abs(eval - median(eval)));
Q2 = mean((pred- eval).^2)  / mean((eval - mean(eval)).^2);
