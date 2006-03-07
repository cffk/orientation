% Convert quaternions to rotations (optionally in turn space).
% $Id$
function x = quattorot(q, turnp),
  if (nargin < 2),
    turnp = false;
  end
  c = q(:,1);
  x = q(:,2:4);
  s = sqrt(sum(x.^2, 2));
  x = x.*repmat((c>=0)*2 - 1, 1, 3);
  c = abs(c);
  a = 2 * atan2(s, c);
  if (turnp),
    a = angletoturn(a);
  end
  x = x.* repmat( a./(s + (s == 0)), 1, 3);
end
