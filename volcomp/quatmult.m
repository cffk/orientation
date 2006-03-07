% multiply two quaternions
% a quaternion is a 4-long row vector
% x and y can be sets of nx and ny quaternions as follows
% nx = ny OR nx = 1 OR ny = 1
% $Id$
function z = quatmult(x,y)
  z = [ x(:,1).*y(:,1)-x(:,2).*y(:,2)-x(:,3).*y(:,3)-x(:,4).*y(:,4), ...
	 x(:,1).*y(:,2)+y(:,1).*x(:,2)+x(:,3).*y(:,4)-x(:,4).*y(:,3), ...
	 x(:,1).*y(:,3)+y(:,1).*x(:,3)+x(:,4).*y(:,2)-x(:,2).*y(:,4), ...
	 x(:,1).*y(:,4)+y(:,1).*x(:,4)+x(:,2).*y(:,3)-x(:,3).*y(:,2)];

