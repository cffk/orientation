% Rotate a vector r by quaternion q
% Either r or q can be sets of (vectors, quaternions) but not both
% $Id$
function s = quattx(q,r)
  s = quatmult(q,quatmult([zeros(size(r,1),1),r],quatconj(q)));
  s = s(:,2:4);
end

