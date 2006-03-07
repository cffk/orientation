% abs of dot product of two quaternions
% a quaternion is a 4-long row vector
% x and y can be sets of nx and ny quaternions as follows
% nx = ny OR nx = 1 OR ny = 1
% $Id$
function z = quatdot(x,y)
  if (size(x,1) == size(y,1)),
    z = abs(sum(x.*y,2));
  elseif (size(x,1) == 1),
    z = abs(sum(repmat(x,size(y,1),1).*y,2));
  else
    z = abs(sum(x.*repmat(y,size(x,1),1),2));
  end
end
