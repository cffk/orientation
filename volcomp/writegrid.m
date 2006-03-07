% Compute and write out a grid file.
% $Id$
function writegrid(delta, num, maxrad, cov),
  [vol rad ind] = dovols(delta,maxrad);
  vol=vol*sum(ind(:,4))/(vol'*ind(:,4));
  f = fopen(['c48u' int2str(num/24) '.grid'], 'wt');
  fprintf(f,'# Orientation set c48u%d, number = %d, radius = %.3f degrees\n',...
	  num/24, num, maxrad);
  fprintf(f,'# $Id$\n');
  fprintf(f,'# For more information, See http://charles.karney.info/orientation/\n');
  fprintf(f, 'format grid\n');
  fprintf(f,'%.6f %.2f %d %d %d %6.3f %8.5f\n',...
	  delta, 0, num, num/24, size(vol,1), maxrad, cov);
  for i=1:size(vol,1),
    fprintf(f, '%2d %2d %2d  %9.7f  %6.3f %3d\n',...
	    ind(i,1), ind(i,2), ind(i,3), vol(i), rad(i), ind(i,4));
  end
  fclose(f);
end
