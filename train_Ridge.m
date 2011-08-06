function w = train_Ridge(X, Y, lambda)
% TRAIN_RIDGE   Traines Ridge regression.
%
%   w = train_Ridge(X, Y, lambda)
%
%   -- input --
%   X:      matrix of input data (number of features x number of examples)
%   Y:      target values
%   lambda: regularisation parameter
%
%   -- output --
%   w:      weight vector of trained Ridge regression
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2007-2011 Georg Zeller, Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2007-2011 Max Planck Society
%


% regularizer
x_dim = size(X, 1);
R = lambda * eye(x_dim, x_dim);
% number of training examples
N = size(X, 2);
Q = X * X';
b = X * Y';
w = (R + Q)\b;
assert(~any(isnan(w)));
