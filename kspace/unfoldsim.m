
close all

%UNFOLD simulation
GAMMA = 26752;  % rad/s/Gauss
debug=0;
tseries=[];

% figure out the constants for the trajectory
dt=0.5e-5;          %   1/(receiver bandwith)
shot_duration = 0.02;   %   total time of the shot.
fov = 64;           % FOV in cm.

% Make the trajectories:
% on the odd images ....
phase = 0;
[traj1 grads]= make_spiral( dt, shot_duration, fov/2, phase);

% do the even images ....
phase = pi;
[traj2 grads]= make_spiral( dt, shot_duration, fov/2, phase);


% Let's create an object to image: 
% object = zeros(64,64);
% object(10:50,10:50)=1;
% object(20:30, 20:30) = 2;
% object=object*10000;

h = read_hdr('sim_0001.hdr');
for n=1:200
    tic
    
    % read each the time points
    str=sprintf('sim_%04d.img',n);
    tmp= read_img2(h,str);
    object = tmp(:,:,1);
    
    %%%%%%
    
    % changing the phase of the spiral allows us to do interleaves
    % reducing the FOV in the spiral design allows us to undersample
    % the object

    if mod(n,2)== 1
        %  kspace1 = ksample(traj1, object, fov);
        %  o1 = ifftshift(ifft2(kspace1))';
        o1 =  ( ksample2(traj1, object, fov) )';      
        
        % save the magnitude images of the original object...
        h.zdim=1;
        str=sprintf('undersamp_%04d.img',n);
        fprintf('\nwriting ... %s',str)
        write_img_data(str,abs(object),h);
        
        str=sprintf('undersamp_%04d.hdr',n);
        write_hdr( str, h );
        
        % save the magnitude images ...
        str=sprintf('alias_sim_%04d.img',n);
        fprintf('\nwriting ... %s',str)
        write_img_data(str,abs(o1),h);
        
        str=sprintf('alias_sim_%04d.hdr',n);
        write_hdr( str, h );
        
        
        % save the phase images too
        str=sprintf('p_alias_sim_%04d.img',n);
        fprintf('\nwriting ... %s',str)
        write_img_data(str,angle(o1)*1000,h);
        
        str=sprintf('p_alias_sim_%04d.hdr',n);
        write_hdr( str, h );
        
        
    else

        %kspace2= ksample(traj2, object,fov);
        %o2 = ifftshift(ifft2(kspace2))';
        o2 =   (ksample2(traj2, object, fov))';
        
        % now that we have a pair, we can combine them, etc.
        %Ktotal = kspace1 + kspace2;
        
        o3 = o1 + o2;
        
        
        % write the magnitude
        str=sprintf('alias_sim_%04d.img',n);
        fprintf('\nwriting ... %s',str)
        write_img_data(str,abs(o2),h);
        
        str=sprintf('alias_sim_%04d.hdr',n);
        write_hdr( str, h );
       
        % write the phase image
        str=sprintf('p_alias_sim_%04d.img',n);
        fprintf('\nwriting ... %s',str)
        write_img_data(str,angle(o2)*1000,h);
        
        str=sprintf('p_alias_sim_%04d.hdr',n);
        write_hdr( str, h );
        
        fprintf('\nsampled k-space for image # %d',n)
        
        % put all the data in a space-time plane:
        o1 = reshape(o1, h.xdim*h.ydim, 1);
        o2 = reshape(o2, h.xdim*h.ydim, 1);
        tseries = [tseries o1 o2];
    end
    toc
    
end

return

%  so now I need to filter the high frequencies ...
fprintf('\nNow applying the filter ...') 
TR = 1;
nyquist = 1/(2*TR);

% the cutoff frequency is expressed in units of Hz.
hcutoff = 0.5 * nyquist;
lcutoff = 0.05 * nyquist;

for n=1:size(object,1)*size(object,2)
    outseries (n,:)= smoothdata2(tseries(n,:), TR, lcutoff,hcutoff,10);
    % this is just for verification:
    if n == 676
        outseries (n,:)= smoothdata2(tseries(n,:), TR, lcutoff,hcutoff, 10,1);
    end
end

% now write the series to file
for n=1:size(outseries,2)
    % magnitude 
    str=sprintf('usim_%04d.img',n);
    fprintf('\n writing .....%s',str)
    write_img_data(str,abs(outseries(:,n)) , h);
    str=sprintf('usim_%04d.hdr',n);
    write_hdr( str, h );
    
    % phase
    str=sprintf('p_usim_%04d.img',n);
    fprintf('\n writing .....%s',str)
    write_img_data(str,1000*angle(outseries(:,n)) , h);
    str=sprintf('p_usim_%04d.hdr',n);
    write_hdr( str, h );  
    
end
%%%%%%%%%%%%%%
%%%  Brad's FFT stuff  ......
% 
% path(path,'/net/fincher/home/bpsutton/Cfiles/')
% path(path,'/net/fincher/home/bpsutton/Matlab/')
%  
% 
% matrix_size = 64;
% 
%  
% %specify kx and ky and time vector, t
%  
% [x,y] = meshgrid([-matrix_size/2:matrix_size/2-1]./matrix_size);
% x = x(:);
% y = y(:);
%  
%  sen = ones(size(x));
%  
% we = zeros(size(x));
% AA1 = mri(real(traj1.')*matrix_size, imag(traj1.')*matrix_size, [0:dt:shot_duration].', we, x, y, sen, 20,64)
% AA2 = mri(real(traj2.')*matrix_size, imag(traj2.')*matrix_size, [0:dt:shot_duration].', we, x, y, sen, 20,64)
% 
% %note last two variables are dummy variables - do not bother changing
%  
% kspace1 = AA1*object(:);
% kspace2 = AA2*object(:);
%  
% img1 = AA1'*kspace1;
% img2 = AA2'*kspace2;
% 
% figure
% subplot(211), imagesc(abs(reshape(img1,64,64)));
% subplot(212), imagesc(abs(reshape(img2,64,64)));
% colormap gray


