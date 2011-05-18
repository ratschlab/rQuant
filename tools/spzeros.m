function A = spzeros(m,n,szmax)
% SPZEROS   Sparse zeros matrix.
%
%   A = spzeros(m,n)
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


if nargin<=1
  A=sparse([],[],[],m,m) ;
  return ;
end 
if nargin<=2
  A=sparse([],[],[],m,n) ;
  return ;
end 
A=sparse([],[],[],m,n,szmax) ;
