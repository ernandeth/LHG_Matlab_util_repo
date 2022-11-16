% spiral3drecon
%
% Luis Hernandez-Garcia @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot
% 3 - it re-grids the k-space data into a cartesian grid
% 4 - reconstructs the image by a 3D FFT
%

dt = 4e-6; %seconds
gamma = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma/2/pi ; % Hz/Gauss

MYDIM = 64;
slthick = 0.8;
myPfile = 'P19456_big_crush'
%myPfile = 'P16384.7'
%myPfile = 'P_nograds'
pfiles = dir('P*.7');
myPfile = pfiles(1).name;

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

%%
Nramp = 500;
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
for Ndel = -160:2:-150
    for Ndel2= 0:5
        figure(1)
        
        % use the k-space trajectory calculated from the gradients,
        % rather than the ideal trajectory
        k2 = gamma*cumsum(g)*dt;
        k3 = k2;
        % spiral out: delays in gradient/acquisition?
        % zero pad first:
        k2 = [zeros(abs(Ndel),3); k2; zeros(abs(Ndel),3)];
        k3 = [zeros(abs(Ndel+Ndel2),3); k3; zeros(abs(Ndel+Ndel2),3)];
        
        % x and y only for the stack of spirals:
        k2(:,1) = circshift(k2(:,1),Ndel,1);
        k2(:,2) = circshift(k2(:,2),Ndel,1);
        % remove the pad
        k2 = k2(abs(Ndel)+1:end-abs(Ndel),:);

        % x and y only for the stack of spirals:
        k3(:,1) = circshift(k3(:,1),Ndel+Ndel2,1);
        k3(:,2) = circshift(k3(:,2),Ndel+Ndel2,1);
        % remove the pad
        k3 = k3(abs(Ndel+Ndel2)+1:end-abs(Ndel+Ndel2),:);
        
        % join the delayed spiral in and out
        k2(end/2+1:end,:) = k3(end/2+1:end,:)
        
        % load k-space trajectory in 2D
        % k2= load('ktraj_cart.txt');
        
        k2 =  k2';
        
        % must resample it to match the data acquisition rate
        % fix this!
        
        ks = zeros(3,ndat);
        fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
        chunk = glen - ndat;
        for n=1:3
            tmp = k2(n,:);
            ks(n,:) = interp1(linspace(0,1,glen), tmp, linspace(0,1,ndat));
        end
        
       
        
        ndat2 = ndat;
        
       

        % if multiple interleaves, concatenate the rotated trajectories
        if nleaves >1
            rotang = 2*pi/nleaves;
            ks_all = [];
            for n=0:nleaves-1
                c = cos(rotang*n); s=sin(rotang*n);
                R = [c s 0; -s c 0; 0 0 1];
                kstmp = R*ks;
                ks_all = [ks_all kstmp]; 
               
            end
            ks = ks_all;
        else
            % throw out ramp
            ks = ks(:, Nramp+1:end-Nramp);
        end
        
        subplot(222)
        plot(ks(1,:), ks(2,:));
        title('Trajectory k-space traj. from grad file ');
        
        %Now do the reconstruction
        % step 1- regrid the k-space data:
        fprintf('\nCreating target cartesian grid ...');
        MYDIM2 = MYDIM*2;
        nslices2 = nslices; %*2;
        
        % Target grid
        Kxymax = max(abs(k2(:)));
        Kxymax = 1* MYDIM/fov;  % Kmax
        Kzmax = 1/(slthick);
        [kx , ky, kz ]= meshgrid( ...
            linspace(-Kxymax,Kxymax,MYDIM2), ...
            linspace(-Kxymax,Kxymax,MYDIM2), ...
            linspace(-Kzmax,Kzmax,nslices2));
        
        % Actually sampled  grid : stack of spirals
        % kx and ky are always the same, but kz increases monotonically
        % this was not recorded in the grad file
        fprintf('\nCreating whole 3D sampled grid ...');
        tmp = repmat(ks,  1, nslices);
        kxs = tmp(1,:);
        kys = tmp(2,:);
        kzs = tmp(3,:);
        ndat3 = ndat2- Nramp*2;
        
        
        kzz = 1*linspace(-Kzmax, Kzmax,nslices);
        for n=1:nslices
            beg = ndat3*(n-1)+1;
            fin = ndat3*n;
            kzs(beg:fin) = kzz(n);
        end
        
       
        
        gkdata = zeros(MYDIM2, MYDIM2, nslices2);
        
        % raw is Nframes x Nslices*ndat x Ncoils
        for n=2 % :nframes
            im_all = 0;
            for c = 1:2:ncoils
                fprintf('\nRegridding frame %d , coil %d...\n', n, c);
                
                tmpk = raw(n,:,c);   % new dims:  Nslices*ndat x 1
                tmpk = reshape(tmpk,ndat*nleaves,nslices);
                
                % throw out ramp
                tmpk = tmpk(Nramp+1:end-Nramp,:);
                
                weights = abs(linspace(0,2,nslices))+1;
                
                for s=1:nslices
                    tmpk(:,s) = tmpk(:,s) * weights(s);
                end
                
                % throw out spiral out:
                %tmpk = tmpk(end/2+1:end,:);
                
                tmpk(isnan(tmpk)) = 0;
                
                
                subplot(221)
                plot(kxs(:)*100+1000)
                hold on
                plot(abs(tmpk(:)));
                hold off
                title('echo train (overlaid on kx wave)')
                drawnow
                
                t_max = find(abs(tmpk(:))==max(abs(tmpk(:))));
                k_center = size(ks,2)*nslices/2 + size(ks,2)/2;
                lag = t_max - k_center
                
                %gkdata = griddata(kxs(:), kys(:), kzs(:), tmpk(:), kx, ky, kz);
                %     gkdata = interp3(kxs(:), kys(:), kzs(:), tmpk(:), kx(:), ky(:), kz(:), 'spline');
                %     F = TriScatteredInterp(kxs(:), kys(:), kzs(:), tmpk(:));
                F = scatteredInterpolant(kxs(:), kys(:), kzs(:), tmpk(:),'linear','none');
                gkdata = F(kx,ky,kz);
                
               
                %gkdata = gkdata(1:2:end, 1:2:end, 1:2:end);
                
                % where the interpolation returned NaN, we must put in zeros
                gkdata(isnan(gkdata)) = 0;
                gkdata = gkdata - mean(gkdata(:));
                subplot(223)
                lightbox(log(abs(gkdata)),[],[]);
                title(sprintf('Log of re-grided k-data, coil %d', c))
                % step 2: do the FFT recon:
                
                im = fft3d(gkdata);
                %im = fftshift(fft2(squeeze(gkdata(:,:,nslices/2+1))));
                
                
                % combine the coils here
                im_all = abs(im_all) + i*abs(im);
                subplot(224)
                lightbox(abs(im_all)); title(sprintf('Reconned, del1=%d, del2=%d', Ndel, Ndel2))
                drawnow
            end
            
            
            
            
        end
        fprintf('\n ... Done! \n');
        %{
        figure(2);
        ov([],abs(im_all),32,32,nslices/2-5,0);
        title(sprintf('Reconned, del1=%d, del2=%d', Ndel, Ndel2))
        drawnow
        
        %}
        imseries(:,:,cnt) = abs(im_all(:,:,nslices/2+1));
        cnt = cnt+1;
    end
end
figure; lightbox(imseries); colormap jet
