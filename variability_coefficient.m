function [Q1 Q2] = variability_coefficient(pred, eval, bg)
% [Q1 Q2] = variability_coefficient(pred, eval)
%
% -- input --
% pred: predicted values
% eval: actual values
% -- output --
% Q1: absolute variability coefficient
% Q2: squared variability coefficient

if 0
Q1 = mean(abs(pred - eval)) / mean(abs(eval - bg));
Q2 = mean((pred - eval).^2) / mean((eval - bg).^2);
else
Q1 = mean(abs(pred - eval)) / mean(abs(eval - median(eval)));
Q2 = mean((pred - eval).^2) / mean((eval - mean(eval)).^2);
end