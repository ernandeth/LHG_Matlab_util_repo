function ktraj=read_ktraj(nfids, doPlots)
% function ktraj=read_ktraj(nfids, doPlots)
%
% nfids is the number of FIDs in each time frame
% doPlots lets you look at each trajectory to make sure everything is OK
%
    g = load('grad.txt');
    M = load('rotmats.txt');
    
    GRADRES = 4e-6; % sampling rate of gradient waveform (s)
    gamma = 267.5*1e2 /(2*pi) ; % Hz/Gauss  
    
    ktraj0 = gamma* cumsum(g,1) * GRADRES;
    
    for n=1:nfids
        m = M( (n-1)*3+1:n*3 , :);
        ktraj{n} = m*ktraj0';
        if doPlots

            subplot(211)
            plot3( ktraj{n}(1,:), ktraj{n}(2,:), ktraj{n}(3,:));
            title('Kspace'); 
            axis([-3 3 -3 3 -3 3])
            axis square
            
            subplot(212)
            plot((m*g')'); title('Gradients');
            
            drawnow; pause(0.25)
        end
    end
    
    
return
    
    