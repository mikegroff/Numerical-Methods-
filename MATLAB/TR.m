% the trapezoidal rule to compute integrals
% fn is the function evaluated at xn
% h = (b-a)/N
function s = TR(fn,h)
N = numel(fn);
s = fn(1) + fn(N) + 2*sum(fn(2:N-1));
s = s*h/2;