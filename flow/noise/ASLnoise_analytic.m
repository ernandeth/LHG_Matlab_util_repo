% this is what the noise input vectors look like
close all; set(gcf,'Position',[10 10 400,600]);
warning off
w=-pi:0.1:pi ;
rho = 0.3;
for noisetype=1:3
   
    switch noisetype
        case 1
            x = ones(size(w));
            %x = w(end:-1:1)/pi;
            title_str='white'
            
        case 2
            x = sqrt(1./w);
            title_str='1/f'
        case 3
            x = ones(size(w));
            x = x ./ (1 - rho * exp(-j.*w));
            x(ceil(length(x)/2))=0;

            title_str='AR(1)'
    end
    subplot(3,1,noisetype)
    plot(w,abs(x))
    title(title_str);
    axis([0 pi 0 3]), axis square

   
    ysimp = x.*(1-exp(-i*w));
    % aliasing: Note that it's very important that you apply the
    % 1-exp(-i*w) factor before you do the aliasing part...
    ysimp = (ysimp + fftshift(ysimp) )/2;   
    ysimp(1:end/4) = 0;
    ysimp(3*end/4: end+1) = 0;
    hold on, 
    plot(w,abs(ysimp),'r')
    
    % running subtraction
    k = (1-exp(-i*(w+pi)))/2;
    yrun = k.*fftshift(x);
    
    plot(w,abs(yrun),'g')
    legend('original','simple','running')

end

