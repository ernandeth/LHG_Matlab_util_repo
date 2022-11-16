function allS = fake_Bfields(N)

allS = [];
mycirc = ones(N);
for n=1:N
    for m=1:N
        if sqrt( (m - N/2-1)^2 + (n - N/2-1)^2) < N/2-2
            mycirc(m,n)=1;
        end
    end
end

% % a test to see what happen to a regioin of extremely low sesitivity
% for n=1:N
%     for m=1:N
%         if sqrt( (m - N/2-1)^2 + (n - N/2-1)^2) < N/8-2
%             mycirc(m,n)=1e-10;
%         end
%     end
% end
% 

for c = 1:8
    t = linspace(0, N, N);
    f = exp(-(t - c* N/8).^2 / (2*pi*N/4)) ;
    
    map = (repmat(f,N,1))' .* mycirc * exp(-i*c*pi/8);
    
    
    %plot(f);
    imagesc(abs(map)); title(sprintf('sensitivity of coil %d', n)),
    drawnow, pause(0.1)
    allS = [allS; map(:)'];
end


return