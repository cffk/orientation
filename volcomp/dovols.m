% Compute volumes and radii of Voronoi cells.
% $Id$
function [vols, rads, inds] = dovols(delta, maxrad),
  turnp = false;
  p = trunccube(delta);
  p0 = cull(p);
  p1 = [];
  s = cubesyms();
  for i = 1:size(s,1),
    p1 = [p1; cull(quatnorm(quatmult(p, s(i,:))), 2 * delta)];
  end
  maxrad = 2 * pi/180*(maxrad+0.02);
  close = cos(maxrad/2);
  vols = [];
  rads = [];
  inds = [];
  for i = 1:size(p0,1),
    p1tx = quatmult(p1, quatconj(p0(i,:)));
    p1tx = p1tx(abs(p1tx(:,1)) >= close, :);
    p1tx = sortrows(p1tx, [-1]);
    p1tx = [1 0 0 0; p1tx(2:size(p1tx,1),:)];
    v = quattorot(p1tx, turnp);
    [vol rad] = voronoivol(v, turnp);
    vols = [vols; vol];
    rads = [rads; rad];
    ind = floor(p0(i,2:4)/p0(i,1)*2/delta+0.5);
    m = 1;
    if (ind(1) ~= 0), m = 2*m; end
    if (ind(2) ~= 0), m = 2*m; end
    if (ind(3) ~= 0), m = 2*m; end
    if (ind(1) == ind(2) & ind(2) == ind(3)),
      m = m;
    elseif (ind(1) > ind(2) & ind(2) > ind(3)),
      m = 6 * m;
    else
      m = 3 * m;
    end
    inds = [inds; [ind m]];
  end
  if (turnp),
    rads = rads * 180/pi * (6*pi)^(1/3);
    vols = vols / (4*pi/3/(24*size(p,1)));
  else
    rads = rads * 180/pi;
    vols = vols / (8*pi^2/(24*size(p,1)));
  end
end
