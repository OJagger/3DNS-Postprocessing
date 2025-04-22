function c = redbluecentered(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end

clim = get(gca, 'Clim');

r = abs(clim(1)/clim(2));
r = r/(r+1);

n1 = floor(r*m);
n2 = floor((1-r)*m);

c1 = redblue(2*n1+1);
c1 = c1(1:n1,:);

c2 = redblue(2*n2+1);
c2 = c2(end-n2+1:end,:);

c = [c1; 1 1 1; c2];

% if (mod(m,2) == 0)
%     % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
%     m1 = m*0.5;
%     r = (0:m1-1)'/max(m1-1,1);
%     g = r;
%     r = [r; ones(m1,1)];
%     g = [g; flipud(g)];
%     b = flipud(r);
% else
%     % From [0 0 1] to [1 1 1] to [1 0 0];
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
% end
% c = [r g b]; 