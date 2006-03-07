% Convert angles to turns.
% $Id$
function t = angletoturn(theta),
  sign = (theta>=0)*2 - 1;
  theta = abs(theta);
  small = theta <= 0.00116863766081;
  t = small .* (0.375750550595608872207534 * theta .* (1 - theta.^2./60)) + ...
      (1 - small) .* ((theta -  sin(theta))/pi).^(1/3);
  t = sign .* t;
end
