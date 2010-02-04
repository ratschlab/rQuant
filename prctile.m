function y = prctile( x, p, dim );

% --- check input
assert( ndims(p) == 2 );
if( size(p,1) == 1 )
  p = p';
end;
assert( size(p,2) == 1 );
assert( all( 0 <= p & p <= 100 ) );
if( ~exist('dim','var') )
  dim = 1;
  if( sum(size(x)>1) == 1 )
    dim = find( size(x) > 1 );
  end;
end;
assert( 1 <= dim & dim <= ndims(x) );

% --- convert percentiles to indices
n = size( x, dim );
j = ( p * n / 100 ) + 0.5;
j0 = floor( j );
nu = j - j0;
assert(all( 0 <= j0 & j0 <= n ));
assert(all( 0 <= nu & nu < 1 ));

% --- special case: 1-dim. data
if( sum(size(x)>1) <= 1 )
  x1 = squeeze( x );
  if( size(x1,1) == 1 )
    x1 = x1';
  end;
  y = prctile1( x1, j0, nu );
  % simulate matlab's dimension inconsistency
  y = y';
  return;
end;

% --- general case: n-dim. data
N = ndims( x );
ydims = size( x );
ydims(dim) = length( p );
y = repmat( nan, ydims );
S = [];
S.type = '()';
S.subs = cell( N, 1 );
% - iterate thru all sub-vectors
% use "subsref" to sequentially single out each vector "x1" from "x"
error( '"prctile" not yet implemented for n-dim. objects "x"' );



% === percentiles of a single vector

function y = prctile1( x1, j0, nu );

assert( ndims(x1) == 2 );
assert( size(x1,2) == 1 );
s = sort( x1 );
s = [ s(1) ; s ; s(end) ];
y = (1-nu) .* s(j0+1) + nu .* s(j0+2);


