GAMMA = 26752  % rad/s/Gauss

%Let's create an object to image:
object = zeros(64,64);
object(16:48,16:48)=1;
subplot(221), imagesc(object), axis('image');

% calculate a field map:
[xx, yy] = meshgrid([-31.5:31.5]);
jj = (xx + i.*yy)';


% calculate magnetic fields from point sources
% abs(B)= 1/r^2 and is shaped like a circle around wire
m1 = cos(angle(jj))./abs(jj);
m2 = cos(angle(jj+8))./abs(jj+8);
m3 = cos(angle(jj+16))./abs(jj+16);

% add the 3 wires
a = m1+m2+m3;
a = m1;



% plot the magnetic field
%subplot(221),colormap(jet);
%mesh(a)

% figure out the color scale and display the field.
amin = min(min(a));
amax = max(max(a));
minmax = [amin,amax]

subplot(222), imagesc(a), axis('image');


% what happens to a voxel in a trajectory?

% figure out the fov stuff
dt=0.5e-5;
samp_time = 0.02;

fov = 64;
kmax = 32/fov; %(samples/ cm)

% spiral gradient waveforms
t=[0:dt:samp_time];
t =t.^(2/3);
spfreq = 2*pi*32/ t(end);


k_spiral=  complex(t .* sin(spfreq*t) , t.*cos(spfreq*t));
k_spiral = k_spiral * kmax/(abs(max(k_spiral))); % normalize

spiral = diff(k_spiral / dt) * (2*pi)/GAMMA;
max_grad = max(real(spiral));

% scale the distorsion field.
a = a* 0.01* max_grad/amax;

% the gradients are the derivative of the field
[gx gy] = gradient(a);
gxy = complex(gx,gy);

% add the signals from all the voxels
signal = zeros(size(k_spiral));
for xi=1:64
    for yi=1:64
                     
        % compute the k-space traj. for each voxel:
        kxy_sus = (GAMMA/(2*pi))* complex(...
            cumsum( gx(xi,yi)* (xi-32) * ones(size(k_spiral))* dt) + ...
            cumsum( gy(xi,yi)* (yi-32) * ones(size(k_spiral))* dt)...
            ) ;
        
        kxy = k_spiral + kxy_sus;
        
        %kxy = k_spiral;
        
        % compute the contribution to the signal:
        signal = signal + ...
            object(xi, yi) * exp(-i * 2*pi .*(real(kxy)*(xi-32) + imag(kxy)*(yi-32)));
        
    end
end

% grid the spiral data and FT
[x,y] = meshgrid([-kmax: kmax/33: kmax]);
kspace=griddata(real(kxy), imag(kxy), signal, x, y); 
kspace(find(isnan(kspace))) = 0;
recon = fftshift(fft2(kspace))';

subplot(223), imagesc(abs(kspace)), title('k space');, axis('image')
subplot(224), imagesc(abs(recon)), title('reconstructed'), axis('image');

% Take a lookk at the trajectories in kspace
xi=1;
while xi>0
    [xi yi] = ginput(1);
    
    xi=round(xi); yi=round(yi);
    [xi yi]
    
    % compute the k-space traj. for each voxel:
    kxy_sus = (GAMMA/(2*pi))* complex(...
            cumsum( gx(xi,yi)* (xi-32) * ones(size(k_spiral))* dt) ,  ...
            cumsum( gy(xi,yi)* (yi-32) * ones(size(k_spiral))* dt)...
            ) ;
    kxy = k_spiral + kxy_sus;
    
    fprintf('\r%d  %d  mag field: %f  gradx: %f  grady: %f\n', ...
        xi,yi, a(xi,yi), gx(xi,yi), gy(xi,yi) );
    
    subplot(224), plot(k_spiral,'r' ), title(sprintf('voxel %d, %d', xi,yi));
    hold on
    subplot(224), plot( kxy_sus,'g' ), title(sprintf('voxel %d, %d', xi,yi));
    subplot(224), plot( kxy ), title(sprintf('voxel %d, %d', xi,yi));
    axis([-1 1 -1 1])
    hold off
end



