% Apply a cubic spread to a function x * (1 + f * x^2)
% $Id$
function y = spread(x, f),
  y = repmat(1+f*sum(x.^2,2),1,3).*x;
end
