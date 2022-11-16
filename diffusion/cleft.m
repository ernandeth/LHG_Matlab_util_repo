function [x, y] = cleft(Nsteps, distance, scale)

% function [x, y] = cleft(Nsteps, distance, scale)
% 
% 

%initialization of variables
y = zeros(size(Nsteps)); 
x = randn(1,Nsteps) * scale;

for n=2:Nsteps

	step = randn(1,1) * scale;

	if (step + y(n-1) >= 0)
		y(n) = step + y(n-1);
	else
		y(n) = y(n-1) + step;
	end
	% note that the particles can not diffuse backwards 

end
tmp = find(y> distance);
if ~isempty(tmp)
	stick = tmp(1);
	y(stick:end) = y(stick);
	x(stick:end) = x(stick);
end

return
