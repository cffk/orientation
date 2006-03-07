% normalized quaternion -- works on a set of quaternions
% $Id$
function y = quatnorm(x)
  y = x./repmat(sqrt(sum(x.^2, 2)), 1, 4);
end



