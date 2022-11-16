option=2;


    Npulses = 24;
    
    myalphas = 10:5:30;
    
    myalphas = 15;
    
if option==1
%% T1 blurring across z direction
    
    data = zeros(1000,2);
    datat = zeros(1000,2);
    
    signal = zeros(Npulses,length(myalphas));
    signalt = zeros(Npulses,length(myalphas));
    
    R1 = 1/1.2;
    dt = 0.001;
    a = 1;
    TR = 30*dt;
    E1 = exp(-TR*R1);
    
    kz = linspace(-12, 12, Npulses);
    if 1
        
        for alpha = myalphas;
            %  initialize an array of pulses
            pulses = linspace(alpha, alpha, Npulses);
            %
            %     for n=2:Npulses
            %         pulses(n) = ...
            %             asin(sin(pulses(1))*tan(pulses(n-1)) / ...
            %             (E1 * sin(pulses(1)) + (1-E1)*tan(pulses(n-1))))
            %     end
            
            
            pulses = [0:Npulses-1].^3 ;
            pulses = alpha + pulses *(90 - alpha)/max(pulses);
            
            %pulses(:) = alpha;
            
            pulses = pulses*pi/180;
            
            % control case:  uninverted spins
            mz=1;
            t=1;
            for p=1:Npulses  % each pulse has a different flip
                
                mx = mz*sin(pulses(p));
                mz = mz*cos(pulses(p));
                
                pnums(p) = round( Npulses/2 + ((-1)^p )*(p/2));
                signal(pnums(p) , a) = mx;
                
                for n=1:30  % compute T1 decay
                    dmz_dt =  R1*(1-mz);
                    mz= mz + dt*dmz_dt;
                    data(t,:) = [mz mx];
                    t = t+1;
                end
            end
            
            % tagged case:  partially inverted spins
            mzt = 1- 2 * exp(-1.2/1.6);
            mzt = 0.95 + 0.05 *mzt;;
            t=1;
            for p=1:Npulses  % each pulse has a different flip
                
                mxt = mzt*sin(pulses(p));
                mzt = mzt*cos(pulses(p));
                
                pnums(p) = round( Npulses/2 + ((-1)^p )*(p/2));
                signalt(pnums(p) , a) = mxt;
                
                for n=1:30  % compute T1 decay
                    dmzt_dt =  R1*(1-mzt);
                    mzt= mzt + dt*dmzt_dt;
                    datat(t,:) = [mzt mxt];
                    t = t+1;
                end
            end
            
            
            a = a+1;
            figure(3) ; plot(data); legend('Mz', 'Mxy') ; title('control');
            figure (4); plot(datat); legend('Mzt', 'Mxt'); title('tag');
            axis([0 800 -1 1])
            
            
            drawnow
        end
        figure(5);
        
        subplot(224);
        plot(kz, signal-signalt);
        axis([-Npulses/2 Npulses/2 0 0.05])
        %legend(num2str(myalphas') )
        ylabel('\DeltaM_{xy}'); xlabel('k_z');
        title('3rd order')
        
        
    end
    
    %
    print -depsc ~/docs/papers/ASL_3D/flip_angles
    figure(7)
    plot(kz,signal,'k');
    hold on
    plot(kz,signalt, 'k--')
    axis([-Npulses/2 Npulses/2 0 0.5])
    ylabel('\DeltaM_{xy}'); xlabel('k_z');
    title('Example: 3rd order schedule starting at 15 deg.')
    legend ('Control', 'Tagged')
    print -depsc ~/docs/papers/ASL_3D/blur_control_tag
    %
    
end

if option==2
    %% Inflow Effects
    % this next simulation shows the blurring resulting from acquiring the data
    % during the uptake.
    width = 2;
    allDelays = 0.3:0.2:1.5
    PSF = zeros(length(allDelays), Npulses);
    delay2 = 0.8;
        
    kz = linspace(-12, 12, Npulses);
    
    d=1;
    for delay1 = allDelays
        PID = 1.5;
        AQtime = 0.5;
        
        DeltaM = zeros(1,Npulses);
        DelataMa = DeltaM;
        
        mydelay1 = delay1;
        
        
        
        for n=1:Npulses
            AQtime = AQtime - 0.02;
            PID = PID + 0.02;
            
            [DeltaM(n) DeltaMa(n)] = ASLconvolution(width, PID, AQtime, mydelay1, delay2)
        end
        
        % tissue signal only:
        tmp = DeltaM-DeltaMa ;
        tmp = tmp/max(tmp);
        % reorder according to k-space acq
        PSF(d, Npulses/2:-1:1) = tmp(1:2:end );
        PSF(d, Npulses/2+1:end) = tmp(2:2:end );
        d=d+1;
    end
    
    figure (5)
    plot(kz, PSF)
    axis([-Npulses/2 Npulses/2 0 1.1])
    title('In-flow Blurring Kernel');
    ylabel('M_{xy} (a. u.)'); xlabel('k_z');
    legend(num2str(allDelays'+delay2) )
    
    print -depsc inflow_blur
    
    save flow_blurring.mat
end
