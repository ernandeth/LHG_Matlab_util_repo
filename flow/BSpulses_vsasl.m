clear
R1a = 1/1.67;
R1g = 1/1.4;
R1w = 1/0.9;
R1csf = 1/2.5;


% in VSASL, it's the time between sat pulse and labeling pulse

tag_delay = 1.25;  % post labeling delay
tag_length = 4.5 - tag_delay - 0.034*(17.5);

AQtime = 1e3*(tag_length + tag_delay);

dt=1e-3;
mcsf(1) = 0;
mg(1) = 0;
mw(1) = 0;
ma(1)= 0;
mat(1)= 0;

SNR = [];

figure(1)
clf


all_BStimes = [];

for BS1_time = [0.05:0.075:1.1]
    for BS2_time =  [0.05:0.075:1.2-BS1_time]
        
        % this is the times at which BS inversion pulses are applied
        % after the labeling time
        BStimes = tag_length + BS1_time + [0 BS2_time];
        BStimes = round(BStimes/dt);
        all_BStimes = [all_BStimes; [BS2_time BS1_time] ];
        
        % Implement Bloch Equation here
        for n=2:AQtime + 1*1e3
            dmg =  (1 - mg(n-1))*R1g;
            mg(n) = mg(n-1) + dmg*dt;
            
            dmw =  (1 - mw(n-1))*R1w;
            mw(n) = mw(n-1) + dmw*dt;
            
            dmcsf =  (1 - mcsf(n-1))*R1csf;
            mcsf(n) = mcsf(n-1) + dmcsf*dt;
            
            dma =  (1 - ma(n-1))*R1a;
            ma(n) = ma(n-1) + dma*dt;
            
            dmat =  (1 - mat(n-1))*R1a;
            mat(n) = mat(n-1) + dmat*dt;
            
            % at the end of the tagging we make sure the arterial spins are
            % tagged
            
            if n == round(tag_length/dt)
                if 1
                    % VSS case : 6850
                    ma(n) = 1;
                    mat(n) = 0;
                    mcsf(n) = mcsf(n)* exp(-28/200);
                    mg(n) = mg(n)* exp(-28/80);;
                    mw(n) = mw(n)* exp(-28/50);;
                    %
                else
                    % VSI case: 6800 and 17268
                    mcsf(n) = -mcsf(n);
                    mg(n) = -mg(n);
                    mw(n) = -mw(n);
                    ma(n) = 1;
                    mat(n) = -1;
                end
            end
            
            % inversion pulse applied to the whole brain and the artery
            for b=1%:length(BStimes)
                tmpBStime = BStimes(b);
                if n == tmpBStime
                    
                    mcsf(n) = -mcsf(n);
                    mg(n) = -mg(n);
                    mw(n) = -mw(n);
                    ma(n) = -ma(n);
                    mat(n) = -mat(n);
                end
            end
            
        end
        %
        
        
        subplot(211)
        
        rectangle('Position', [AQtime -1 600  2],'FaceColor',[0.5 0.5 0.5] )
        hold on
        plot(mg,'k'); grid on
        plot(mw,'y'); grid on
        plot(mcsf,'g');
        %plot(ma,'b')
        %plot(mat,'r')
        %plot(abs(ma - mat),'c')
        legend('Gray Matter', 'White Matter', 'CSF','blood','tagged blood','difference','Location','NorthWest')
        hold off
        
        Mtotal1 = ma + 50*mw + 50*mg + 50*mcsf;
        Mtotal2 = mat + 50*mw + 50*mg + 50*mcsf;
    
        
        subplot(212)
        rectangle('Position', [AQtime -1 600  2],'FaceColor',[0.5 0.5 0.5] )
        hold on
        plot(abs(Mtotal1) - abs(Mtotal2))
        %}
        hold off
        %rectangle('Position', [AQtime -1 600  2],'FaceColor',[0.5 0.5 0.5] )
        
        tmp = (abs(ma) - abs(mat)) ./ (mw + mg + mcsf);
        
        tmp = mcsf; % look at the CSF alone
        tmp = mean(tmp(AQtime:AQtime+10));
        SNR = [SNR tmp];
        
        drawnow
        
    end
end

figure(2)
plot((SNR))
