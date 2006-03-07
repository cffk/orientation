% compute the volume and radius of a voronoi cell corresponding to the
% first point.  In flat space we need a = b = 0.  However rotation and
% turn space are not flat.  To compensat, a spreads out the initial
% points and then b spreads out the vertices of the Voroinoi cells.  a
% and b were determined by minimizing the error relatve to an "exact"
% quadrature for delta > 0.05.  Estimated max error is (delta/12)^2
% and the error is worst for the (odd-shaped) boundary cells.  The
% error for the interior cells is much less.
% $Id$
function [volume, rad] = voronoivol(x, turnp),
  if (nargin < 2),
    turnp = false;
  end
  if (turnp)
    f=3.55;
    a = f*1/20;
    b = (1-f)*2.6/20;
  else
    % Adjustment for conversion to rotation space.  (6*pi)^(1/3) is the
    % linear scaling.  1/60 is the coefficient of the cubic term in
    % the Taylor expansion of the turn formula.
    % a = a/(6*pi)^(2/3)-1/60;
    % b = b/(6*pi)^(2/3);
    % Here are optimized values based on c48u8649 (close to previous ones)
    f = -0.8442;
    a = f*(-0.00961);
    b = (1-f)*(-0.025);
  end
  [v, c] = voronoin(spread(x,a));
  [k, volume] = convhulln(spread(v(c{1},:),b));
  rad = max(sqrt(sum(v(c{1},:).^2,2)));
end
