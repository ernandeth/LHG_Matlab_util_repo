% this version is intended to simulate conditions at 7T
%
% Batch file for pseudoCASL simulations.  COntains a bunch of cases
% execute a chunk of code at a time


doBatch=1;
sigs=[];
phis = [];
close all;
%
array_Gz_ss = linspace(0.1,2.5,15);
array_Gz_ave = linspace(0,0.15, 15);

sigs1 = zeros(length(array_Gz_ave), length(array_Gz_ss));
n=1;
figure(1)
for Gz_ss = array_Gz_ss
    m=1;
    for Gz_ave = array_Gz_ave
        
        tag_loc = 0;
        flip_ang = 22.5;
        isTag = 1;
        vel = 35;
        slomo =0;
        eta = 0.4;   % works for all locs if eta = 0.25!!
        t_ramp = 0.02;  % ramp time in ms.
        t_rfrepeat = 1.2;
        
        %Gz_ss = 1;  % slice select gradient in G/cm.
        % Gz_ave = 0.05;  % this is the average gradient between pulses
        Gz_err = 0;
        pCASL09
        
        tagSignal = (M(end,3));
        sigs1(m,n) = tagSignal;
        
        m = m+1;
        
    end
    
    n = n+1;
end

figure(4)
subplot(3,2,1)
imagesc(sigs1)
colorbar
ylabel('Mean Gz');
xlabel('Max. Gz');

set(gca,'Xtick',[1,5,10,15]);
set(gca,'Ytick',[1,5,10,15]);

set(gca,'YtickLabel',(array_Gz_ave([1,5,10,15])));
set(gca,'XtickLabel',(array_Gz_ss([1,5,10,15])));
%
%%


array_vel = linspace(10,100,15);
array_flips =linspace(5,65,15);
sigs2 = zeros(length(array_vel), length(array_flips));
n=1;
figure(1)
for m = 1:length(array_vel)
    for n = 1:length(array_flips)
        
        isTag = 1;
        tag_loc = 0;
        
        flip_ang = array_flips(n);
        vel = array_vel(m);
        
        fprintf('\nvelocity: %02f , flip: %02f', vel, flip_ang);
        slomo =0;
        eta = 0.4;   % works for all locs if eta = 0.25!!
        t_ramp = 0.02;  % ramp time in ms.
        
        Gz_ss = 1;  % slice select gradient in G/cm.
        Gz_ave = 0.05;  % this is the average gradient between pulses
        Gz_err = 0;
        pCASL09
        
        tagSignal = (M(end,3));
        sigs2(m,n) = tagSignal;
    end
    
end


figure(4)
subplot(3,2,2)
imagesc(sigs2)
colorbar
ylabel('Velocity');
xlabel('Flip Angle');

set(gca,'Xtick',[1,5,10,15]);
set(gca,'Ytick',[1,5,10,15]);

set(gca,'XtickLabel',(array_flips([1,5,10,15])));
set(gca,'YtickLabel',(array_vel([1,5,10,15])));
%%

array_vel = linspace(10,100,15);
array_Gz_ave = linspace(0,0.15,15);

sigs3 = zeros(length(array_vel), length(array_Gz_ave));
n=1;
figure(1)
for m = 1:length(array_vel)
    for n = 1:length(array_Gz_ave)
        
        isTag = 1;
        tag_loc = 0;
        
        flip_ang = 22.5;
        vel = array_vel(m);
        
        slomo =0;
        eta = 0.4;   % works for all locs if eta = 0.25!!
        t_ramp = 0.02;  % ramp time in ms.
        
        Gz_ss = 1;  % slice select gradient in G/cm.
        Gz_ave = array_Gz_ave(n);  % this is the average gradient between pulses
        Gz_err = 0;
        pCASL09
        
        tagSignal = (M(end,3));
        sigs3(m,n) = tagSignal;
    end
    
end


figure(4)
subplot(3,2,3)
imagesc(sigs3)
colorbar
ylabel('Velocity');
xlabel('Mean Gz');

set(gca,'Xtick',[1,5,10,15]);
set(gca,'Ytick',[1,5,10,15]);

set(gca,'XtickLabel',(array_Gz_ave([1,5,10,15])));
set(gca,'YtickLabel',(array_vel([1,5,10,15])));

%%

array_flips =linspace(5,65,15);
array_trep = linspace(0.6,1,15);

sigs4 = zeros(length(array_flips), length(array_trep));
n=1;
figure(1)
for m = 1:length(array_trep)
    for n = 1:length(array_flips)
        
        isTag = 1;
        tag_loc = 0;
        
        
        vel = 35;
        
        slomo =0;
        eta = 0.4;   % works for all locs if eta = 0.25!!
        t_ramp = 0.02;  % ramp time in ms.
        
        t_rfrepeat = array_trep(m);
        flip_ang = array_flips(n);
        
        Gz_ss = 1;  % slice select gradient in G/cm.
        Gz_ave = 0.05;  % this is the average gradient between pulses
        Gz_err = 0;
        pCASL09
        
        tagSignal = (M(end,3));
        sigs4(m,n) = tagSignal;
    end
    
end


figure(4)
subplot(3,2,4)
imagesc(sigs4)
colorbar
ylabel('Repeat Time');
xlabel('Flip Angle');

set(gca,'Xtick',[1,5,10,15]);
set(gca,'Ytick',[1,5,10,15]);

set(gca,'XtickLabel',(array_flips([1,5,10,15])));
set(gca,'YtickLabel',(array_trep([1,5,10,15])));

%%
%%

array_flips =linspace(5,65,15);
array_Gz_ave = linspace(0,0.15, 15);
sigs5 = zeros(length(array_flips), length(array_Gz_ave));
n=1;
figure(1)
for m = 1:length(array_Gz_ave)
    for n = 1:length(array_flips)
        
        isTag = 1;
        tag_loc = 0;
                
        vel = 35;
        
        slomo =0;
        eta = 0.4;   % works for all locs if eta = 0.25!!
        t_ramp = 0.02;  % ramp time in ms.
        
        t_rfrepeat = 1.2; 
        
        flip_ang = array_flips(n);
                
        Gz_ss = 1;  % slice select gradient in G/cm.
        Gz_ave = array_Gz_ave(m);  % this is the average gradient between pulses
        Gz_err = 0;
        pCASL09
        
        tagSignal = (M(end,3));
        sigs5(m,n) = tagSignal;
    end
    
end


figure(4)
subplot(3,2,5)
imagesc(sigs5)
colorbar
ylabel('Gz mean');
xlabel('Flip Angle');

set(gca,'Xtick',[1,5,10,15]);
set(gca,'Ytick',[1,5,10,15]);

set(gca,'XtickLabel',(array_flips([1,5,10,15])));
set(gca,'YtickLabel',(array_Gz_ave([1,5,10,15])));

return

%
% parms that work:
%     flip_ang = 22.5;
%     isTag = 1;
%     vel = 40;
%     slomo =0;
%     eta = 0.25;   % works for all locs if eta = 0.25
