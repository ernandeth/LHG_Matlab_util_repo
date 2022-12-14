% abadd.m - add gain in absolute units
%
% Usage: h = abadd(type, style, th, g, rays, width)
%
% type = 1,2 for polar, azimuthal
% style = line style of added gain curve, e.g., '--r'
%
% rest of arguments as in ABPOL and ABZ
% 
% use h to add legends: h1 = abz(phi, g1);
%                       h2 = abadd(2, '--', phi, g2);
%                       legend([h1,h2], 'gain 1', 'gain 2');

% S. J. Orfanidis - 1997 - www.ece.rutgers.edu/~orfanidi/ewa

function h = abadd(type, style, th, g, rays, width)

if nargin==0, help abadd; return; end
if nargin<5, rays = 30; end
if nargin<6, width = 1; end
   
if type == 1,
   x = g .* sin(th);                             % polar plot
   y = g .* cos(th);
   h = plot(x, y, style, 'LineWidth', width);
   hold on;
   plot(-x, y, style, 'LineWidth', width);
else
   x = g .* cos(th);                             % azimuthal plot
   y = g .* sin(th);
   h = plot(x, y, style, 'LineWidth', width);
   hold on;
   plot(x, -y, style, 'LineWidth', width);
end

