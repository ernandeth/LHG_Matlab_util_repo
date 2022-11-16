clear all
close all

Num_Sections = 4; %number of sections in laminar flow profile
Ra= 0.5; % Radius of Artery in cm
V_max = 60; %(cm/s)

off_start= 0;
off_step = 10; %0.02;
off_end = 500;

grad_err_start = -0.3;
grad_err_step = 0.02;
grad_err_end = 0.3;


Y=off_start:off_step:off_end;
X=grad_err_start:grad_err_step:grad_err_end;


As=[]; %area of each disk (...( ( ( (A1)A2)A3)A4)... ) => cm^2
Vs=[];
A_temp = 0;
for n = Num_Sections:-1:1
    ar_disk = pi*(Ra/n)^2;
    ar_slice = ar_disk - A_temp;
    As=[As ar_slice];
    A_temp = ar_disk;
    
         
    Vs=[Vs V_max*(1 - (Ra/n - Ra/Num_Sections)/Ra )^2];

end

for ii=1:Num_Sections
    C(ii)=As(ii)*Vs(ii)/(sum(As.*Vs));
end


IEs=[];
cntr_flip=1;

for off_res_peak = off_start:off_step:off_end
    
    cntr_eta=1;
    for Gz_err_amp=grad_err_start:grad_err_step:grad_err_end

        IE_temp=[];
        for vel=Vs


            showplot = 0;
            pcolor='red';
            %hold on;
            %%%%%%%%%%%%% TAGGING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tag_loc = 0;
            %vel = 60;
            %%%%%%%%%%%%% TAGGING PULSE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flip_ang = 30; % in Degrees
            eta = 0.3; % (A+ - A-)/ A+
            isTag = 1;

            flip_ang = flip_ang * pi /180;
            RF_spacing = 0.80 ; % time between RF pulses in ms.
            t_ramp = 0.02;  % ramp time in ms.
            Gz_ss = 0.6;  % slice select gradient in G/cm.
            t_Gz_ss = 0.5 ; % ms
            t_Gz_rephaser = 0.200 ; % ms

            xtra_phase2 =0; %(rad)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%% GRADIENT ERROR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
            Gradient_error = 1; % 0: No gradient error | 1 : gradient error
            Gradient_error_shape = 2; %1: global  2: local
            %Gz_err_amp = 0.1;     % amplitude of gradient error
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%OFF-RESONANCE ERROR PARAMETERS%%%%%%%%%%%%%%%%
            off_res_error = 1;
            %off_res_peak = 350; % in Hz
            off_res_width = 3; %cm - width of off-resonance pattern
            off_shape=1; % 1:rect , 2:hanning - shape of off-resonance pattern
            error_loc = 3 ; % 1:before; 2:first 1/4; 3:middle; 4:first 3/4; 5: last
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fixerror = 0; % 0: No fixing| 1: eta update| 2: prephasing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            slomo = 0;


            pCASL_offres_sim2
            IE_temp= [IE_temp (1-M(end,3))/2*exp(250/1664)];
        end


        IEs (cntr_flip,cntr_eta)=  sum(IE_temp.*C);
        cntr_eta=cntr_eta +1;
    end
    cntr_flip = cntr_flip +1;
end


figure;
[C,h]=contour(X,Y,IEs);
ylabel('off resonance (HZ)');
xlabel('Gradient error (G/cm)');
set(h,'ShowText','on');
title(sprintf('OffPeak/W=%2.2d(HZ)/%2.2g(cm)-MaxVel=%2.2g(cm/s)-fixmethod=%2.2d',off_res_peak,off_res_width,V_max,fixerror));



figure;
[C,h]=contour(X,Y,IEs,[0.8 0.7]);
ylabel('off resonance (HZ)');
xlabel('Gradient error (G/cm)');
set(h,'ShowText','on');
title(sprintf('OffPeak/W=%2.2d(HZ)/%2.2g(cm)-MaxVel=%2.2g(cm/s)-fixmethod=%2.2d',off_res_peak,off_res_width,V_max,fixerror));


figure;
mesh(X,Y,IEs);
ylabel('off resonance (HZ)');
xlabel('Gradient error (G/cm)');
zlabel('inversion efficiency');
colorbar;
title(sprintf('OffPeak/W=%2.2d(HZ)/%2.2g(cm)-MaxVel=%2.2g(cm/s)-fixmethod=%2.2d',off_res_peak,off_res_width,V_max,fixerror));

figure;
imagesc(X,Y,IEs);
ylabel('off resonance (HZ)');
xlabel('Gradient error (G/cm)');
zlabel('inversion efficiency');
colorbar;
title(sprintf('OffPeak/W=%2.2d(HZ)/%2.2g(cm)-MaxVel=%2.2g(cm/s)-fixmethod=%2.2d',off_res_peak,off_res_width,V_max,fixerror));

figure;
[C,h]=contour(X,Y,IEs);
ylabel('off resonance (HZ)');
xlabel('Gradient error (G/cm)');
colorbar;
title(sprintf('OffPeak/W=%2.2d(HZ)/%2.2g(cm)-MaxVel=%2.2g(cm/s)-fixmethod=%2.2d',off_res_peak,off_res_width,V_max,fixerror));

savename= sprintf('GRAD-OFF-OffPeakW%d-%g-MaxVel%g-fixmethod%d',off_res_peak,off_res_width,V_max,fixerror);


str=['save ' savename];
eval(str);

% plot([eta_start:eta_step:eta_end],IEs,'-*');
% grid on;
% title(sprintf('Vmax=%2.2g (cm/s) , Flip Angel=%2.2g',V_max,flip_ang2));
% ylabel('Inversion Efficiency');
% xlabel('\eta ');