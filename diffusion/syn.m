function net = syn(Nsteps,Nspins,scale,grad, ratio, dt)

% function net = syn(Nsteps,Nspins,scale,grad, ratio,dt)
% 
% returns the NET xy-magnetization vector of a population of spins composed
% of two subpopulations : restricted diffusion and inrestricted
%
% Nsteps:  numberof steps in the spin history
% Nspins: total number of spins in the ensemble
% scale:  scaling factor that determines the size of the diffusion jumps
% gamma:  gyromagnetic constant
% grad:   a vector with a magnetic field gradient waveform
% ratio: the ratio of restricted/unrestricted spins
% dt: the duration of the individual diffusion jumps

GAMMA = 4258*2*pi;  % rad/s / Gauss
% vesicle spins:
    radius = 4e-8;     %source:  Brad's pharmacology book
    count=0;
    M=[]; N=[]; 
    Mx=[]; My=[];
    Nx=[]; Ny=[];
    for jc = 1:round(Nspins*ratio)
        fprintf('\rspin number %d', jc)
        % compute the trajectory of the particles
        [xx yy] = bounce(Nsteps, radius, scale);
        
        % compute the phase gain from the gradient
        phi = GAMMA* dt* cumsum(grad.*xx);
        % compute the Mag. vectors
        magnet = ones(size(phi)).* exp((i* phi));
        M = [M; magnet]; 
        count=count+1;
        Mx = [Mx ; xx];
        My= [My; yy];
    end
    


    % cleft spins:
    distance = 2e-7;
   
    for jc = round(Nspins*ratio)+1:Nspins
        fprintf('\rspin number %d', jc)
        % compute the trajectory of the particles
        [xx yy] = cleft(Nsteps, distance, scale);
        
        % compute the phase gain from the gradient
        phi = GAMMA* dt* cumsum(grad.*xx);
        % compute the Mag. vectors
        magnet = ones(size(phi)).* exp((i* phi));
        
        N = [N; magnet];
        count=count+1;
        Nx = [Nx; xx];
        Ny= [Ny; yy];
    end
    
    MM = sum(M,1);
    NN = sum(N,1);
    
    if isempty(NN)
        NN=zeros(size(MM));
    end
    
    if isempty(MM)
        MM=zeros(size(NN));
    end
    net = (MM+NN);%.*exp(-t/T2);
    fprintf('\n')
    
    str = sprintf('save vecs%.2f.mat Mx My Nx Ny MM NN', ratio);
    eval(str);
    
 return
