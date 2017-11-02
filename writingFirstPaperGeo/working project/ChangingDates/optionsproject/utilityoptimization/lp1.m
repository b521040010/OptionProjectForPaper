function ret = lp1( x )
%LP1 (Log of |x| plus 1 )*sgn(x). Allows a number from -inf to inf to be
%    shown on a log-like scale
ret = log( abs(x)+1) .* sign(x);
end

