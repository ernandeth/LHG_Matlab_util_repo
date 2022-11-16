function [allS ] = read_BFields(basestr, N)
% this is the function that reads in the .mat files with the simulated B
% fields

allS = [];
msk = zeros(N,N);
dx = 128/N;  % resampling
for r=1:N
    for c=1:N
        if sqrt((r-N/2)^2  + (c-N/2)^2)< N/2-20
            msk(r,c)=1;
        end
    end
end


for n=1:8
    str = [basestr num2str(n) '.mat']
    
    load(str);
    Bx = abs(Bx); 
    By = abs(By); 
    B = complex(Bx,By);
    % normalize the sensitivities
    maxB = max(abs(B(:)));
    B = B/maxB;
    B = B(1:dx:end, 1:dx:end);    
    B = B .*msk;
    % Note  .. this doesn't make a difference on the G map. 
    % Now, try masking out the weirdness on the edge:

    figure(1); imagesc(abs(B)); drawnow,    title(sprintf('sensitivity of coil %d', n)),  pause(0.1);
    allS = [allS; B(:)'];
end
return