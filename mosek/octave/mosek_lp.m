function [x,obj,INFO,lambda] = mosek_lp (x0, q, A, b, lb, ub, A_lb, A_in, A_ub, num_threads)
% [x,obj,INFO,lambda] = mosek_lp (x0, q, A, b, lb, ub, A_lb, A_in, A_ub, num_threads)
%
% Written (W) 2008-2009 Fabio De Bona, Gunnar Raetsch
% Copyright (C) 2008-2009 Max Planck Society
%
% This file is part of the Mosek Wrapper for Octave.
%
% -*- texinfo -*-
% @deftypefn {Function File} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} lp (@var{x0}, @var{H}, @var{q}, @var{A}, @var{b}, @var{lb}, @var{ub}, @var{A_lb}, @var{A_in}, @var{A_ub})
% Solve the quadratic program
% @iftex
% @tex
% $$
%  \min_x {1 \over 2} x^T H x + x^T q
% $$
% @end tex
% @end iftex
% @ifnottex
%
% @example
%      min 0.5 x'*H*x + x'*q
%       x
% @end example
%
% @end ifnottex
% subject to
% @iftex
% @tex
% $$
%  Ax = b \qquad lb \leq x \leq ub \qquad A_{lb} \leq A_{in} \leq A_{ub}
% $$
% @end tex
% @end iftex
% @ifnottex
%
% @example
%      A*x = b
%      lb <= x <= ub
%      A_lb <= A_in*x <= A_ub
% @end example
% @end ifnottex
%
% @noindent
% using a null-space active-set method.
%
% Any bound (@var{A}, @var{b}, @var{lb}, @var{ub}, @var{A_lb},
% @var{A_ub}) may be set to the empty matrix (@code{[]}) if not
% present.  If the initial guess is feasible the algorithm is faster.
%
% The value @var{info} is a structure with the following fields:
% @table @code
% @item solveiter
% The number of iterations required to find the solution.
% @item info
% An integer indicating the status of the solution, as follows:
% @table @asis
% @item 0
% The problem is feasible and convex.  Global solution found.
% @item 1
% The problem is not convex.  Local solution found.
% @item 2
% The problem is not convex and unbounded.
% @item 3
% Maximum number of iterations reached.
% @item 6
% The problem is infeasible.
% @end table
% @end table
% @end deftypefn

   if ~(nargin == 4 || nargin == 6 || nargin == 9 || nargin==10)
      print_usage ();
   else

      % Checking the initial guess (if empty it is resized to the
      % right dimension and filled with 0)
      if (isempty (x0))
         x0 = zeros (n, 1);
      elseif (length (x0) ~= n)
         error ('lp: the initial guess has incorrect length');
      endif

      % Linear penalty.
      if (length (q) ~= n)
         error ('lp: the linear term has incorrect length');

      endif

      % Equality constraint matrices
      if (isempty (A) || isempty(b))
         n_eq = 0;
         A = zeros (n_eq, n);
         b = zeros (n_eq, 1);
      else
         [n_eq, n1] = size (A);

         if (n1 ~= n)
            error ('lp: equality constraint matrix has incorrect column dimension');
         endif

         if (length (b) ~= n_eq)
            error ('lp: equality constraint matrix and vector have inconsistent dimension');
         endif

      endif

      % Bound constraints
      Ain = zeros (0, n);
      bin = zeros (0, 1);
      n_in = 0;
      if (nargin > 4)
         if (~ isempty (lb))
            if (length(lb) ~= n)
               error ('lp: lower bound has incorrect length');
            else
               Ain = [Ain; eye(n)];
               bin = [bin; lb];
            endif
         endif

         if (~ isempty (ub))
            if (length (ub) ~= n)
               error ('lp: upper bound has incorrect length');
            else
               Ain = [Ain; -eye(n)];
               bin = [bin; -ub];
            endif
         endif
      endif

      % Inequality constraints
      if (nargin > 6)
         [dimA_in, n1] = size (A_in);
         if (n1 ~= n)
            disp (' dimension are '), disp(n1), disp(n)
            error ('lp: inequality constraint matrix has incorrect column dimension');
         else
            if (~ isempty (A_lb))
               if (length (A_lb) ~= dimA_in)
                  error ('lp: inequality constraint matrix and lower bound vector inconsistent');
               else
                  Ain = [Ain; A_in];
                  bin = [bin; A_lb];
               endif
            endif

            if (~ isempty (A_ub))
               if (length (A_ub) ~= dimA_in)
                  error ('lp: inequality constraint matrix and upper bound vector inconsistent');
               else
                  Ain = [Ain; -A_in];
                  bin = [bin; -A_ub];
               endif
            endif
         endif
      endif

      P = zeros(0,0)
   
      if ( nargin > 9 )
         [x] = __mosek_qp__ (x0, P, q, A, b, lb, ub, A_lb, A_in, A_ub, num_threads);
      elseif ( nargin > 6 )
         [x] = __mosek_qp__ (x0, P, q, A, b, lb, ub, A_lb, A_in, A_ub);
      elseif ( nargin > 4 )
         [x] = __mosek_qp__ (x0, P, q, A, b, lb, ub);
      else
         [x] = __mosek_qp__ (x0, P, q, A, b);
      endif

   endif

endfunction
