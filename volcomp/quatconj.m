% conjugate of a quaternion -- works on a set of quaternions
% $Id$
function y = quatconj(x)
  y = [x(:,1), -x(:,2:4)];
end
