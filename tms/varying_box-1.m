% straight.m								
%
% This program models the electric field generated when electrodes
% are placed at two sides of a square strip of conducting paper. The paper is 
% modeled by a grid. 
% The student can set the size of the box and the tolerance
voltage = input('What is the voltage of the nongrounded side?');
boxside = input('Number of units on side of box? ');

% Set the maximum allowable change in value for a point between iterations
tolerance = input('tol ');	

% The next line initializes the grid, setting all points to zero volts
old = zeros(boxside,boxside);
old(boxside,:) = voltage + zeros(1,boxside);  % set one side of box to 9 volts
old(boxside,:) = [0:boxside-1]*voltage/boxside;

equilibrium = 0;                % a variable used to check for equilibrium
icount = 0;                     % an iteration counter
disp('Done with initialization of boundary conditions.');

new = old;

ih1 =         figure('Units','normalized','Position',[0 .5 .45 .5]);
ih2 =         figure('Units','normalized','Position',[.5 .5 .45 .5],'KeyPressFcn','squareone');

while (equilibrium < 1)
    for m = 2:boxside-1	
        for n = 1:boxside
            new(m,n) = (old(m-1,n)+old(m+1,n))/2;
        end
    end
    icount = icount + 1;		% add one to the iteration counter
    figure(ih1);
    
    colormap(hot);
    pcolor(new);
    shading flat;
    hold on;
    contour(new,8,'k');
    hold off;
    title('Straight');
    colorbar
    drawnow;    
    
    % Position the figure in the upper right hand corner of the screen
    figure(ih2);
    surf(new);
    colormap(hot);
    colorbar;
    drawnow;    
    if all(new-old<=tolerance)	% test for equilibrium
        equilibrium = 1;
        disp('number of iterations');
        disp(icount);
        
        close;		% Clean Up previous figures
        close;
        % Position the figure in the upper right hand corner of the screen
        figure(ih1);
        
        colormap(hot);
        pcolor(new);
        shading flat;
        hold on;
        contour(new,8,'k');
        hold off;
        title('Straight');
        colorbar
        drawnow;    
        
        % Position the figure in the upper right hand corner of the screen
        figure(ih2);
        surf(new);
        colormap(hot);
        colorbar;
        drawnow;    
    end
    old = new;
end	