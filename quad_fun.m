function f = quad_fun(x, S1, S2, S3)
% f = quad_fun(x, S1, S2, S3)
%
% -- input --
% x: variable (scalar)
% S1: coefficient of quadratic term
% S2: coefficient of linear term
% S3: constant
%
% -- output --
% f: function value


f = x^2*S1 + x*S2 + S3;
