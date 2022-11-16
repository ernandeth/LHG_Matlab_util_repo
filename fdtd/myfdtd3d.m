%***********************************************************************
%     3-D FDTD code with PEC boundaries
%***********************************************************************
%     modified by Luis Hernandez to include variable maps of conductivity
%     and permeability
%     date march 2005
%
%     Program author: Susan C. Hagness
%                     Department of Electrical and Computer Engineering
%                     University of Wisconsin-Madison
%                     1415 Engineering Drive
%                     Madison, WI 53706-1691
%                     608-265-5739
%                     hagness@engr.wisc.edu
%
%     Date of this version:  February 2000
%
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a three-Nensional
%     Cartesian space lattice comprised of uniform cubic grid cells.
%     
%     To illustrate the algorithm, an air-filled rectangular cavity 
%     resonator is modeled.  The length, width, and height of the 
%     cavity are 10.0 cm (x-direction), 4.8 cm (y-direction), and 
%     2.0 cm (z-direction), respectively.
%
%     The computational domain is truncated using PEC boundary 
%     conditions:
%          ex(i,j,k)=0 on the j=1, j=jb, k=1, and k=kb planes
%          ey(i,j,k)=0 on the i=1, i=ib, k=1, and k=kb planes
%          ez(i,j,k)=0 on the i=1, i=ib, j=1, and j=jb planes
%     These PEC boundaries form the outer lossless walls of the cavity.
%
%     The cavity is excited by an additive current source oriented
%     along the z-direction.  The source waveform is a differentiated 
%     Gaussian pulse given by 
%          J(t)=-J0*(t-t0)*exp(-(t-t0)^2/tau^2), 
%     where tau=50 ps.  The FWHM spectral bandwidth of this zero-dc-
%     content pulse is approximately 7 GHz. The grid resolution 
%     (dx = 2 mm) was chosen to provide at least 10 samples per 
%     wavelength up through 15 GHz.
%
%     To execute this M-file, type "fdtd3D" at the MATLAB prompt.
%     This M-file displays the FDTD-computed Ez fields at every other
%     time step, and records those frames in a movie matrix, M, which 
%     is played at the end of the simulation using the "movie" command.
%
%***********************************************************************

 clear

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

%***********************************************************************
%     Grid parameters
%***********************************************************************

N=30;       %number of grid cells in x-direction
M=N-1;
   

is=N/2;       %x location of z-directed current source
js=N/2;       %y location of z-directed current source

kobs=12;      % observation plane

dx=0.002;          %space increment of cubic lattice
dt=dx/(2.0*cc);    %time step

nmax=100;          %total number of time steps

%***********************************************************************
%     Differentiated Gaussian pulse excitation
%***********************************************************************

rtau=150.0e-12;
tau=rtau/dt;
ndelay=3*tau;
%ndelay=0;%3*tau;
srcconst=-dt*3.0e+11;

%***********************************************************************
%     Material parameters
%***********************************************************************
%
%eps=1.0;
%sig=0.0;        
%***********************************************************************
%     Material parameter maps
%***********************************************************************
epsz=ones(N,N,N)*epsz;
sig=zeros(N,N,N);
%sig(10:15,10:15,10:15)=1e-10;
mu = muz * ones(N,N,N);

%***********************************************************************
%     Updating coefficients
%***********************************************************************
% 
% ca=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
% cb=(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));
% da=1.0;
% db=dt/muz/dx;
%***********************************************************************
%     Updating coefficients
%***********************************************************************

Ca = 1.0-(dt.*sig)./(2.0*epsz);
Cb = dt/(epsz*2*dx);

Da = 1.0 - dt .* sig./(mu);
Db = dt./(mu*2*dx);

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(N,N,N);
ey=zeros(N,N,N);
ez=zeros(N,N,N);
hx=zeros(N,N,N);
hy=zeros(N,N,N);
hz=zeros(N,N,N);


%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************
src= srcconst.*([1:500]-ndelay).*exp(-(([1:500]-ndelay).^2./tau^2));
subplot(223), plot(src);



for n=1:nmax
      
%***********************************************************************
%     Update electric fields
%***********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %          Let's try a Cartesian Grid
%            I only update elements from 2:N-1  - note M=N-1 to abbreviate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
% Ex = Ca*Ex + Cb( dHy/dz - dHz/dy); 
ex(1:N, 2:M, 2:M)= Ca(1:N, 2:M, 2:M) .* ex(1:N, 2:M, 2:M) - ...
                   Cb (1:N, 2:M, 2:M).* ...
                   (hy(1:N, 2:M, 3:M+1) - hy(1:N, 2:M, 1:M-1) - ...
                    hz(1:N, 3:M+1, 2:M) + hz(1:N, 1:M-1, 2:M));  %- Jx(2:M, 2:M, 2:M) ...
                               
% Ey = Ca*Ey + Cb( dHx/dz - dHz/dx); 
ey(2:M, 1:N, 2:M) = Ca(2:M, 1:N, 2:M) .* ey(2:M, 1:N, 2:M) - ...
                    Cb(2:M, 1:N, 2:M) .* ...
                    (hx(2:M, 1:N, 3:M+1) - hx(2:M, 1:N, 1:M-1) - ...
                     hz(3:M+1, 1:N, 2:M) + hz(1:M-1, 1:N, 2:M)) ;  %- Jy(2:N,1:N,2:N) ...
                     
% Ez = Ca*Ey + Cb( dHy/dx - dHx/dy); 
ez(2:M, 2:M, 1:N) = Ca(2:M, 2:M, 1:N) .* ez(2:M, 2:M, 1:N) - ...
                    Cb(2:M, 2:M, 1:N) .* ...
                    (hy(3:M+1, 2:M, 1:N) - hy(1:M-1, 2:M, 1:N) - ...
                     hx(2:M, 3:M+1, 1:N) + hx(2:M, 1:M-1, 1:N)) ; %- Jz(2:N,1:N,2:N) ...
                 
% this will go away when I uncomment the J terms:
ez(is:is+1, js:js+1, :) = ez(is:is+1, js:js+1,:) + ...
                0.0001*srcconst.*(n-ndelay)*exp(-((n-ndelay)^2/tau^2));

%***********************************************************************
%     Update magnetic fields
%***********************************************************************

% Hx = Da*Hx + Db*(dEy/dz + dEz/dy - Mx*dt ) 
hx(1:N, 2:M, 2:M) =  Da(1:N, 2:M, 2:M) .* hx(1:N, 2:M, 2:M) + ...
                     Db(1:N, 2:M, 2:M) .* ...
                     (ey(1:N, 2:M, 3:M+1) - ey(1:N, 2:M, 1:M-1) - ...
                      ez(1:N, 3:M+1, 2:M) + ez(1:N, 1:M-1, 2:M)); 
                                        % - Mx*dt/mu);
                    
% Hy = Da*Hy + Db*(dEz/dx - dEx/dz - My*dt) 
hy(2:M, 1:N, 2:M) = Da(2:M, 1:N, 2:M) .* hy(2:M, 1:N, 2:M) + ...
                    Db(2:M, 1:N, 2:M) .* ...
                    (ez(2:M, 1:N, 3:M+1) - ez(2:M, 1:N, 1:M-1) - ...
                     ex(3:M+1, 1:N, 2:M) + ex(1:M-1, 1:N, 2:M)) ;
                                        % - My*dt);
                    
% Hz = Da*Hz + Db*(dEx/dy - dEy/dx - Mz*dt) 
hz(2:M, 2:M, 1:N) = Da(2:M, 2:M, 1:N) .* hz(2:M, 2:M, 1:N) + ...
                    Db(2:M, 2:M, 1:N) .* ...
                    (ex(2:M, 3:M+1, 1:N) - ex(2:M, 1:M-1, 1:N) - ...
                     ey(3:M+1, 2:M, 1:N) + ey(1:M-1, 2:M, 1:N)) ;
                     % - Mz*dt);
                    

%***********************************************************************
%     Visualize fields
%***********************************************************************

if mod(n,2)==0;

    timestep=int2str(n);
    tview(:,:)=ez(:,:,N/2);
    sview(:,:)=ez(:,N/2,:);
    subplot(221),%;
  %   plot(squeeze(ex(:, N/2+1,N/2)))
    imagesc(tview')
%     shading flat;
    caxis([-1.0 1.0]);
    colorbar;
    axis image;
    title(['Ez( z=N/2), time step = ',timestep]);
    
    subplot(222),
%     plot(squeeze(ex(N/2,:,N/2)))
%     title(['time step = ', timestep])
    imagesc(sview');
%     shading flat;
    caxis([-1.0 1.0]);
    colorbar;
    axis image;
    title(['Ez(y=N/2), time step = ',timestep]);
    
    drawnow

    subplot(224), 
    plot(squeeze(ez(:,N/2,  N/2)))
    xlabel('x position')
    %axis([0 N -100 100])
    nn=n/2;
    %M(:,nn)=getframe(gcf,rect);

end;

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

%movie(gcf,M,0,10,rect);
