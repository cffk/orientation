% Produce a body centered cubic mesh with indices n.  E.g., bcc(-5:5).
% $Id$
function p = bcc(n),
  [x y z] = meshgrid(n, n, n);
  m = (bitand(abs(x),1) == bitand(abs(y),1)) & ...
      (bitand(abs(y),1) == bitand(abs(z), 1));
  p = [z(m) x(m) y(m)];
end

