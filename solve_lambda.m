function lambda = solve_lambda(R1, R2, R0)
% lambda = solve_lambda(R1, R2, R0)
%
% -- input --
% R1: coefficient of cubic term
% R2: coefficient of constant term
% R0: coefficient of quartic term
%
% -- output --
% lambda: positive real solution of the quartic function


% solves 2*R0*lambda^4 + R1*lambda^3 - 2*R2= 0
if nargin==3
  % solution from http://www.wolframalpha.com/input/?i=solve%282+R0+lambda^4%2BR1+lambda^3-2+R2%29
  tmp_lambda(1) = -R1/(8*R0)+sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))/2-sqrt(R1^2/(8*R0^2)+(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))-(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0)-R1^3/(32*R0^3*sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))))/2;
  tmp_lambda(2) = -R1/(8*R0)+sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))/2+sqrt(R1^2/(8*R0^2)+(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))-(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0)-R1^3/(32*R0^3*sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))))/2;
  tmp_lambda(3) = -R1/(8*R0)-sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))/2-sqrt(R1^2/(8*R0^2)+(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))-(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0)+R1^3/(32*R0^3*sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))))/2;
  tmp_lambda(4) = R1/(8*R0)-sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))/2+sqrt(R1^2/(8*R0^2)+(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))-(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0)+R1^3/(32*R0^3*sqrt(R1^2/(16*R0^2)-(8*R2)/(3^(1/3)*(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3))+(-9*R1^2*R2+sqrt(3)*sqrt(27*R1^4*R2^2+4096*R0^3*R2^3))^(1/3)/(2*3^(2/3)*R0))))/2;

  fidx = find(imag(tmp_lambda)==0 & tmp_lambda>0);
  try 
    assert(length(fidx)==1);
  catch
    error('choice of regularisation parameters not appropriate');
  end
  lambda = tmp_lambda(fidx);
else
  lambda = (2^(1/3)*R2^(1/3))/R1^(1/3);
end