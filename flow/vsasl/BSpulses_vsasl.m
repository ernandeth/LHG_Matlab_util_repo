clear
close all

R1a = 1/1.65;
R1=1/1.4;
R1csf = 1/3;

dt=1e-3;

TR = 4.5
RO_time = 0.6;
t_delay = 1.4;

t_adjust = TR - RO_time-t_delay;


do7t = 0;

if do7t
    R1a = 1/2.3;
    R1=1/2;
    R1csf = 1/4.3;
    t_delay = 0.3;
    t_adjust = 2.5;
    BS2_time = [0.005:0.02:t_delay];
    BS2_time = 0.5
    tag_end = 1;
end


%BS1_time= linspace(0.1, t_delay-0.1, 8) ;
BS1_time = 0.1;

BS2_time = 0.25;

AQtime = 1e3*(t_adjust + t_delay);

% asl is the observed signal for each combination
asl=zeros(8);
bkgnd=asl;

for n1 = 1:length(BS1_time)
    
    %BS2_time = linspace(0.1, t_delay - BS1_time(n1)-0.1, 8) ;
    
    for n2 =1:length(BS2_time)
        
        
        b = t_adjust + BS1_time(n1);
        c = b + BS2_time(n2);
        b = floor(b/dt);
        c = floor(c/dt);
        
        %%
        mcsf(1) = 0;
        mcsft(1) = 0;
        m(1) = 0;
        ma(1)= 0;
        mat(1)= 0;
        
        for n=2:AQtime
            dm =  (1 - m(n-1))*R1;
            m(n) = m(n-1) + dm*dt;
            
            dmcsf =  (1 - mcsf(n-1))*R1csf;
            mcsf(n) = mcsf(n-1) + dmcsf*dt;
            
            dmcsft =  (1 - mcsft(n-1))*R1csf;
            mcsft(n) = mcsft(n-1) + dmcsft*dt;
            
            dma =  (1 - ma(n-1))*R1a;
            ma(n) = ma(n-1) + dma*dt;
            
            dmat =  (1 - mat(n-1))*R1a;
            mat(n) = mat(n-1) + dmat*dt;
            
            
            
            
            %labeling pulse
            if n == round(t_adjust/dt)
                mcsf(n) = -mcsf(n);
                mcsft(n) = -mcsft(n);
                m(n) = -0.9*m(n);
                ma(n) = 0;
                mat(n) = -0.9*mat(n);
            end
            
            % BGS inversion pulse 1 applied to the whole brain and the artery
            if n == b
                mcsf(n) = -mcsf(n);
                mcsft(n) = -mcsft(n);
                m(n) = -0.8*m(n);
                ma(n) = -0.8*ma(n);
                mat(n) = -0.8*mat(n);
            end
            
            % BGS inversion pulse 2 applied to the whole brain and the artery
            if n == c
                mcsf(n) = -mcsf(n);
                mcsft(n) = -mcsft(n);
                m(n) = -m(n);
                ma(n) = -ma(n);
                mat(n) = -mat(n);
            end
            
        end
        
        figure(1)
        rectangle('Position', [(t_adjust)/dt, -1 , 0.01/dt , 2],'FaceColor',[0.5 0.5 0.5] )
        
        
        
        %         plot(mcsf,'g')
        %         plot(mcsft,'y');
        %         plot(ma,'b')
        %         plot(mat,'r')
        plot(m,'b'); grid on ;
        hold on;
        plot(mat-ma,'c')
        plot(mcsf,'k')
        legend('tissue','asl con1rast','csf con1rast')
        %        legend('tissue','CSF','CSF tag','blood','tagged blood','difference','Location','NorthWest')
        
        
        asl(n1,n2) =  mat(end-2)-ma(end-2);
        bkgnd(n1,n2) =  m(end-2);
        bkgnd2(n1,n2) =  mcsf(end-2) ;
        drawnow; pause(0.2);
        hold off
    end
end

  %
% figure
% subplot(211)
% imagesc(asl); title('ASL signal'); colorbar
% ylabel('Transit Time')
% xlabel('in1erpulse BS2_time')
%
% subplot(212)
% imagesc(abs(bkgnd)); title('Bcknd Signal'); colorbar
% ylabel('Transit Time')
% xlabel('in1erpulse BS2_time')
% legend('ASL signal','tissue')
