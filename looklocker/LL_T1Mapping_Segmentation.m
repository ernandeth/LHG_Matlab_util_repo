%  This script is for estimating the T1 map and
%  the fraction maps for the GM, WM and CSF using the LL images
clear all
close all
 
 
Filename = 'vol_e4164_09_16_116_0029.nii'; 

Gen_File = strtok(Filename, '.');

TR = 350 %  ms
% Flip angle is 10 degrees
alf = 10*pi/180;

% Nominal T1 values
T1_W = 700
T1_G = 1500
T1_C = 3500


% Water fraction in brain regions
WF_W = 0.73;
WF_G = 0.89;
WF_C = 1;
 


 vars = load_nii(Filename);
 nframes =  vars.hdr.dime.dim(5);
 nslice = vars.hdr.dime.dim(4);


 IRim = double(vars.img);
  
 
 F_map_W = zeros(vars.hdr.dime.dim(2), vars.hdr.dime.dim(3), vars.hdr.dime.dim(4));
 F_map_G = zeros(vars.hdr.dime.dim(2), vars.hdr.dime.dim(3), vars.hdr.dime.dim(4));
 F_map_C = zeros(vars.hdr.dime.dim(2), vars.hdr.dime.dim(3), vars.hdr.dime.dim(4));
 T1_map = zeros(vars.hdr.dime.dim(2), vars.hdr.dime.dim(3), vars.hdr.dime.dim(4));
 M0_map = zeros(vars.hdr.dime.dim(2), vars.hdr.dime.dim(3), vars.hdr.dime.dim(4));
 
 for x = 1:vars.hdr.dime.dim(2)
   
   h = waitbar(x/vars.hdr.dime.dim(2), 'Please Wait....');
    close(h)

    for y = 1:vars.hdr.dime.dim(3)
        
        for z = 1:vars.hdr.dime.dim(4)
                
        if (IRim(x,y,z,1) > 500) 
        tcurve = (squeeze(IRim(x,y,z,:)))';

        % Finding the zero crossing point in  the LL signal
        gg = find(tcurve == min(tcurve));
        tcurve(1:gg-1) = -1*tcurve(1:gg-1);
        tcurve(gg:end) = tcurve(gg:end);

        % Finding M0 and T1
        LBA = 0;
        T1 = [100,tcurve(1)];
        
        % Calculating the sampling time points 
        ti = zeros(1,nframes);
        for i = 1:nframes
            ti(i) = 17.5+TR*i;
        end

        % Fitting the data points to the equation
        [T1, fval1, exit_flag]=fminsearchbnd(@(T1_EST) LL_T1(T1_EST,tcurve, ti,alf,TR),T1);

        
        T1_1 = T1(1);
        M0 = T1(2);
        
        T1_map(x,y,z)= T1_1;
        M0_map(x,y,z) = M0;
        
      
        
        W = zeros(nframes,1);
        G = zeros(nframes,1);
        C = zeros(nframes,1);
        S = zeros(nframes,1);
        

        
        T1_W_m = T1;
        T1_W_m(1) = T1_W;
        W = LL_T1_test(T1_W_m,ti,alf,TR);
        
        T1_G_m = T1;
        T1_G_m(1) = T1_G;        
        G = LL_T1_test(T1_G_m,ti,alf,TR);
        
        T1_C_m = T1;
        T1_C_m(1) = T1_C;               
        C = LL_T1_test(T1_C_m,ti,alf,TR);
        
        S = tcurve';
   
        
        
        AMAT = [W',G',C'];
        
        AMAT_N = AMAT'*AMAT;
        S_N = AMAT'*S;
        
        AMAT_N_I = inv(AMAT_N);
        alfM = AMAT_N_I*S_N;
        

                
        
        F_map_W(x,y,z) = alfM(1);
        F_map_G(x,y,z) = alfM(2);
        F_map_C(x,y,z) = alfM(3);
       

        end
       end

    end
 end
 
 
 
 vars.hdr.dime.dim(1) = 3;
 vars.hdr.dime.dim(5) = 1;
 vars.original.hdr = vars.hdr;
 vars.hdr.dime.bitpix = 64;
 vars.hdr.dime.datatype = 64;
 
 
 vars.img = double(T1_map);
 maxim = max(max(max(T1_map)));
 minim = min(min(min(T1_map)));
 vars.hdr.dime.glmax = maxim;
 vars.hdr.dime.glmin = minim;
 save_nii(vars, strcat('T1_map_',Gen_File,'.nii'));
 
 
 vars.img = double(F_map_W);
 maxim = max(max(max(F_map_W)));
 minim = min(min(min(F_map_W)));
 vars.hdr.dime.glmax = maxim;
 vars.hdr.dime.glmin = minim;
 save_nii(vars, strcat('F_map_W_',Gen_File,'.nii'));
 
 vars.img = double(F_map_G);
 maxim = max(max(max(F_map_G)));
 minim = min(min(min(F_map_G)));
 vars.hdr.dime.glmax = maxim;
 vars.hdr.dime.glmin = minim;
 save_nii(vars, strcat('F_map_G_',Gen_File,'.nii'));
 
  
 vars.img = double(F_map_C);
 maxim = max(max(max(F_map_C)));
 minim = min(min(min(F_map_C)));
 vars.hdr.dime.glmax = maxim;
 vars.hdr.dime.glmin = minim;
 save_nii(vars, strcat('F_map_C_',Gen_File,'.nii'));
 

 
 