% spiral3drecon
%
% Luis Hernandez-Garcia @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot
% 3 - reconstructs the image by Jeff's NUFFT
%

dt = 4e-6; %seconds
gamma_rad = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma_rad/2/pi ; % Hz/Gauss

MYDIM = 64;


pfiles = dir('P*.7');
myPfile = pfiles(1).name;

doRots = 0;
doGoldenAngle = 0;
doTraj= 0;
delta_x = 7; % plane offset for Trajectory measurement
delta_y = 7; % plane offset for Trajectory measurement

fprintf('\nBegin CP NUFFT recon of spiral trajectory for %s', myPfile)
%run('~hernan/matlab/irt/setup')

% load the raw data:
% Nframes x Nslices*ndat x Ncoils
[raw scaninfo ] = read_raw_3d(myPfile ,0 );

nslices = scaninfo.nslices;
nechoes = scaninfo.npr;
nframes = scaninfo.nphases;
ndat = scaninfo.ndat;
ncoils = scaninfo.ncoils;
nfids = nechoes*nframes*nslices;
fov = scaninfo.opfov;
MYDIM = scaninfo.opxres;
nleaves = scaninfo.npr;
slthick = scaninfo.slthick/10;

fov3 = [fov fov nslices*slthick];
dim3 = [MYDIM MYDIM nslices];

figure(1)
subplot(221)
plot(squeeze((abs(raw(1,:,:)))));
title('Echo train: Frame 1, all coils');
drawnow



% Now reshape into a 3D matrix
% kdata = reshape(kdata,  ndat, nslices, nframes);

fprintf('\nLoading grads...');
g = load('grad.txt');
glen = size(g,1);
aqwindow = dt*glen;
dt2 = aqwindow/ndat;
%%
Nramp =0;
slord = [];
m=1
for n=1:nslices/2
    slord(m) =  nslices/2 + (n);
    slord(m+1) =  nslices/2 +1-(n);
    m=m+2;
end
slord
imseries = [];
cnt = 1;
for Ndel = 0 %42
   
    for Ndel2= 0
        figure(1)
        
        % load k-space trajectory in 2D
        % k2= load('ktraj_cart.txt');
        
        % use the k-space trajectory calculated from the gradients,
        % rather than the ideal trajectory
        k2 = gamma*cumsum(g)*dt;  % this is in cm-1 units
        k3 = k2;
        
        % including delays in gradient/acquisition?
        % spiralin and spiralout may have different delays
        % do the shifts separately, then rejoin.
        
        
        % spiralin:
        % zero pad first:
        k2 = [zeros(abs(Ndel),3); k2; zeros(abs(Ndel),3)];
        
        % shift the x and y only for the stack of spirals:
        k2(:,1) = circshift(k2(:,1),Ndel,1);
        k2(:,2) = circshift(k2(:,2),Ndel,1);
        % remove the pad
        k2 = k2(abs(Ndel)+1:end-abs(Ndel),:);
        
        % spiralout:
        % zero pad first:
        k3 = [zeros(abs(Ndel+Ndel2),3); k3; zeros(abs(Ndel+Ndel2),3)];
        % shift the x and y only for the stack of spirals:
        k3(:,1) = circshift(k3(:,1),Ndel+Ndel2,1);
        k3(:,2) = circshift(k3(:,2),Ndel+Ndel2,1);
        % remove the pad
        k3 = k3(abs(Ndel+Ndel2)+1:end-abs(Ndel+Ndel2),:);
        
        % join the delayed spiral in and out
        k2(end/2+1:end,:) = k3(end/2+1:end,:)
        %}
        % tranpose to make things work
        k2 =  k2';
                
        % adjust for measured GRADIENT ERROR ???
        % k2 = k2 *0.68;

        % The error may be here:
        %
        % must resample the trajectory to match the data acquisition rate
        ks = zeros(3,ndat);
        fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
        chunk = glen - ndat;
        for n=1:3
            tmp = k2(n,:);
            %ks(n,:) = interp1(linspace(0,1,glen), tmp, linspace(0,1,ndat));
            ks(n,:) = tmp(1:ndat);
        end
        
        
        % use a thinner slice, and
        % collect a reference frame WITHOUT gradients
        
        
        
        % if multiple interleaves,
        % generate the rotated version and then concatenate the rotated trajectories
        if nleaves >1
            rotang = pi/nleaves;
            ks_all = [];
            for n=0:nleaves-1
                c = cos(rotang*n); s=sin(rotang*n);
                R = [c -s 0; s c 0; 0 0 1];
                kstmp = R*ks;
                kstmp = kstmp(:, Nramp+1:end-Nramp) ;
                ks_all = [ks_all kstmp];
            end
            ks = ks_all;
            %Nramp = 0;
            
        end
        
        % calculate kx, ky kz locations for all the stack of spirals
        % kx and ky are always the same, but kz increases monotonically
        % this was not recorded in the grad file
        fprintf('\nCreating whole 3D sampled grid ...');
        
        % first throw out ramp
        if nleaves < 2
            ks = ks(:, Nramp+1:end-Nramp) ;
        end
        ndat3 = ndat - 2*Nramp;
        
        subplot(222)
        plot(ks(1,:), ks(2,:));
        title('Trajectory k-space traj. from grad file ');
        
        if doRots==0
            % copy the xy trajectory to all platters
            tmp = repmat(ks,  1, nslices);
            kxs = tmp(1,:);
            kys = tmp(2,:);
            kzs = tmp(3,:);
            
            % calculate the kz components
            Kzmax = 1/(slthick);
            kzz = 1*linspace(-Kzmax, Kzmax,nslices);
            for n=1:nslices
                beg = nleaves*ndat3*(n-1)+1;
                fin = nleaves*ndat3*n;
                kzs(beg:fin) = kzz(n);
                % kzs(beg:fin) = kzz(slord(n));
            end
            % rejoin the components into a single matrix
            ks = [kxs' kys' kzs'];
        else
            tmpks = [];
            xrot_angles = pi*(0:nslices-1)/nslices;
            if (doGoldenAngle)
                 xrot_angles = 137.5*pi/180*[0:nslices-1];
            end
            
            for s= xrot_angles
                xi = s *pi/(nslices);
                cxi = cos(xi); sxi = sin(xi);
              
                
%                 Rx =[cxi 0   sxi ;
%                      0   1   0 ;
%                    -sxi  0  cxi];

                  Rx =[ 1  0   0 ;
                        0 cxi  sxi ;
                        0 -sxi cxi];
                
                tmp =  Rx*ks;
                plot3(tmp(1,:), tmp(2,:), tmp(3,:));
                axis(Kxymax*[-1 1 -1 1 -1 1]);
                drawnow
                
                tmpks = [tmpks tmp];
                
            end
            ks = tmpks';
            
        end
        
        
        plot3(ks(:,1), ks(:,2), ks(:,3));
        
        
        MYDIM2 = MYDIM*2;
        
        Kxymax = max(abs(k2'));
        Kxymax = max(Kxymax(1:2));
        %Kxymax = MYDIM/fov;  % Kmax
        fprintf('\nCreating target cartesian grid ...');
        %MYDIM2 = MYDIM*2;
        
        % Target grid
        %Kxymax = max(abs(k2(:)));
        Kxymax = 1* MYDIM/fov/2;  % Kmax
        Kzmax = 1/(slthick);
        [kx , ky, kz ]= meshgrid( ...
            linspace(-Kxymax,Kxymax,MYDIM), ...
            linspace(-Kxymax,Kxymax,MYDIM), ...
            linspace(-Kzmax,Kzmax,nslices));

       
        
        % including a filter to un-do the T2 effect
        T2 = 0.150;
        
      
        %}
        
        % raw is Nframes x Nslices*ndat x Ncoils
        for f=1 %:nframes 
            im_all = 0;
            im = zeros(MYDIM,MYDIM,nslices);
            for c = 1:ncoils
                fprintf('\nReconning frame %d , coil %d...\n', n, c);
                
                tmpk = raw(f,:,c);   % new dims:  Nslices*ndat*nleaves x 1
                tmpk = reshape(tmpk,ndat*nleaves,nslices);
                b = reshape(tmpk, ndat, nleaves, nslices);
                b = b(Nramp+1:end-Nramp,:,:);
                tmpk = reshape(b,(ndat-2*Nramp)*nleaves, nslices);
                % filter to account for T2 effects along Z axis
                %weights = 2*abs(linspace(0,2,nslices))+1;
                %weights = exp(linspace(0,nslices*55, nslices)/200);
                
                % throw out spiral out?
                % tmpk = tmpk(end/2+1:end,:);
                
                tmpk(isnan(tmpk)) = 0;
                
                figure(1)
        
                subplot(221)
                plot(ks(:,2)/Kxymax)
                hold on
                plot(abs(tmpk(:))/max(abs(tmpk(:))));
                hold off
                title('echo train (overlaid on kx wave)')
                drawnow
                
                % gridding recon
                for s=1:nslices
                    F = scatteredInterpolant(ks(1:end/2,1), ks(1:end/2,2), tmpk(:,s),'linear','none');
                    gkdata = F(kx(:,:,s),ky(:,:,s));
                    gkdata(isnan(gkdata)) = 0;
                    lightbox(log(abs(gkdata)),[],[]);
                    im(:,:,s) = fftshift(fft2(fftshift(gkdata)));
                end
                lightbox(log(abs(gkdata)),[],[]); title('Log k-space image (re-grided)')
                
                %im = nufft2d(ks(1:end/2,1:2), tmpk(:,2), 0.01*fov3(1:2), dim3(1:2));
                
                %xcp = nufft3d(ks, tmpk(:), fov3, dim3);  
                figure(2)
                
                % coil combination:
                im_all = abs(im_all) + i*abs(im);
                %}
                figure
                lightbox(abs(im_all)); title(sprintf('Reconned, del1=%d, del2=%d', Ndel, Ndel2))
                if doTraj
                    figure(1)
                    
                    buf = raw(f,:);
                    buf = reshape(buf,ndat,nslices);
                    buf = buf(:,nslices/2);
                    %buf = buf(:,3);
                    
                    if n==1 % this is a regular frame
                        whole = buf;
                        
                    end
                    
                    if n==2 % this is the reference frame: no grads.
                        ref = buf;
                        subplot(222)
                        hold off
                        plot(phase(whole)/gamma_rad/dt2/1e3)
                        hold on
                        plot(unwrap(angle(ref)))
                        hold off
                        title('frames 1 and 2: signal phases')
                        legend('both grads on', 'no grads (ref)')
                    end
                    
                    if n==3 % x gradient only
                        buf = buf ./ref;
                        kx_measured = unwrap(angle(buf))/delta_x/2/pi;
                        subplot(223)
                        hold off
                        plot(kx_measured)
                        hold on
                        plot(ks(1:ndat,1));
                        hold off                        
                        legend('kx measured',  'kx nominal')
                        title(sprintf('Kx (frame 3)'))

                    end
                    if n==4
                        buf = buf ./ref;
                        ky_measured = unwrap(angle(buf))/delta_y/2/pi;
                        subplot(224)
                        hold off
                        plot(ky_measured)  
                        hold on
                        plot(ks(1:ndat, 2));
                        legend(  'ky measured', 'ky nominal')
                        title(sprintf('Ky (frame 4)'))
                        hold off
                        
                    end
                end
            end
            
            
            fprintf('\n ... Done! \n');
            %
            figure(2);
            ov([],abs(im_all), MYDIM/2, MYDIM/2, nslices/2, 0);
            title(sprintf('CP Reconned, del1=%d, del2=%d', Ndel, Ndel2))
            drawnow
            %}
        end
        %
        %
        %imseries(:,:,cnt) = abs(im_all(:,:,nslices/2+1));
        %cnt = cnt+1;
        %}
    end
end
%%

figure; lightbox(abs(imseries)); colormap jet
%{
 figure(2);
 zfov = nslices*slthick/2;
 for n=1:nslices
%       subplot(221)
%      set(gca,'Position',[ 0.05    0.55    0.5  0.5*fov/zfov ]);
%      subplot(222)
%      set(gca,'Position',[ 0.55    0.55    0.5  0.5*fov/zfov ]);
%      subplot(223)
%      set(gca,'Position',[ 0.05    0.05    0.5  0.5 ]);
     ov([],abs(im_all),MYDIM/2,MYDIM/2,n,0);
%    drawnow
     pause(0.5)
 end
%}