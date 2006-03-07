% Return a set of bcc points in trucated cube
% $Id$
function v = trunccube(delta),
  n = ceil((sqrt(2) - 1)/(delta/2));
  v = bcc(-n:n);
  v = v * delta/2;
  m = (abs(v(:,1)) < sqrt(2) - 1) & ...
      (abs(v(:,2)) < sqrt(2) - 1) & ...
      (abs(v(:,3)) < sqrt(2) - 1) & ...
      (abs(v(:,1)) + abs(v(:,2)) + abs(v(:,3)) < 1);
  v = v(m,:);
  v = quatnorm([repmat(1, size(v,1),  1), v]);
end
