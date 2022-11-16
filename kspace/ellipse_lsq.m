function result = ellipse_lsq(parms, xdata,  ydata)

x0 = parms(1);
y0 = parms(2);
a = parms(3);
b = parms(4);

doFig=1;

if nargin > 2
	% if we are fitting the data with an ellipse, we use the 
	% x coordinates provided by the users
	ydata = reshape(ydata,length(ydata),1);
    xdata = reshape(xdata, length(xdata),1);
    
	% compute top part of ellipse
	y = y0 + abs(b*sqrt( 1 - ((xdata - x0)/a).^2));
	% compute bottom half of elllipse
	yneg = y0 - abs(b*sqrt( 1 - ((xdata - x0)/a).^2));
	
	
	% calculate cost function.  
	% We use the top or bottom, depending which is 
	% smaller
	cost = zeros(size(ydata));
	for count=1:length(ydata)
		if (ydata(count) - yneg(count)) < (ydata(count) - y(count))
			cost(count) = ydata(count) -  y(count);
		else
			cost(count) = ydata(count) -  yneg(count);
		end
	end
	result = cost;
			
    if doFig
        plot(xdata,ydata,'*')
        hold on
        plot(xdata,y)
        plot(xdata,yneg)
        hold off
        drawnow
        pause
    end
else
	% if we are not fitting a curve, generate the x-coordinates
	% given the a parameter and the ellipse location
    % essentially ignore the input for x
	if isempty(xdata)
		x = linspace(x0-a, x0+a,20);
	end
	% compute top part of ellipse
	y = y0 + b*sqrt( 1 - ((x - x0)/a).^2);
	% compute bottom half of elllipse
	yneg = y0 - b*sqrt( 1 - ((x - x0)/a).^2);
	% put together:
	y = [y ,yneg(end:-1:1)];
	x = [x x(end:-1:1)];
	result = [x' y'];
end

return
