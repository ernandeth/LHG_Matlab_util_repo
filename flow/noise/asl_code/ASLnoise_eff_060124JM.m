% function ASLnoise_eff_112905JM
% script to illustrate the properties of the different subtraction methods
% in ASL - turboCASL in event related experiments.
%
% uses simulated data with AR(1)+WN noise
%
% 5 differencing matrices.
% calculates power, efficiency, MSE, bias of the beta estimate as a
% function of SNR
%
%
% % set up the options here:
% DesignType = 2;  % ASL with no BOLD:  1=blocked, 2=ER Fixed, 3=ER randomized
%                  % ASL with BOLD:-1=blocked,-2=ER Fixed,-3=ER randomized
% ResultsFile = 'Blocksim_results';
% doSpectra = 0;
% showPlots = 0;
% warning off
% rho = 0.9;
% VarAR = 0.1;  % Ratio of White Noise variance to AR variance (0=AR only)
% VarWN = 10
% TR = 1.4;
% nyq = 1/(2*TR);
% tlen = 258;
% NITER = 1;
% doNewMatrix = 1;  % generate a new design matrix or read it from a file?

% % determine a range of noise levels - SNRs..
% 
% nlevels = [0.1:0.1: 2];  Set nlevels to the values of SNR that you want
% Want SNR's from 0.1 thorugh 2
% 

disp(['Working on : ' ResultsFile])


noBOLD = DesignType>0; %this indicates if a BOLD regressor is used or not

% Make the differencing matrices
D1 = eye(tlen);
D2 = zeros(tlen/2,tlen);
for count=1:tlen/2
    D2(count,count*2-1)=1;
    D2(count,count*2)=-1;
end
D3 = zeros(tlen-1,tlen);
for count=1:tlen-1
    D3(count,count)=(-1)^(count-1);
    D3(count,count+1)=(-1)^(count);
end
D4 = zeros(tlen-2,tlen);
for count=1:tlen-2
    D4(count,count)=(-1)^(count-1);
    D4(count,count+1)=2*(-1)^(count);
    D4(count,count+2)=(-1)^(count-1);
end
% make up the sinc-interpolation matrix.
klen = 41;
t = linspace(-round(klen/2),round(klen/2),klen);
skernel = (sinc(t-0.5));
sk = zeros(1,length(skernel)*2);
sk(2:2:end) = skernel;
D5 = zeros(tlen + 2*klen);
for row=2:2:tlen+klen
    D5(row,row-1:row-2 + 2*klen ) = sk;

end
for row=3:2:tlen+klen
    D5(row,row-1:row-2 + 2*klen ) = -sk;
end

D5 = D5(klen:end-klen-1, 2*klen:end-klen+1);
for row=1:tlen
    % normalize the kernel
    D5(row,:) = D5(row,:) / abs(sum(D5(row,:)));
    % stick in the ones to represent the non-interpolated samples
    D5(row,row) = -(-1)^row;
end
D5 = D5(1:tlen, 1:tlen);

nDtype = 6; % D1-D5 OLS plus D1 GLS
Spower_SNR = zeros(nDtype,length(nlevels));
eff_SNR = zeros(nDtype,length(nlevels));
bias_SNR = zeros(nDtype,length(nlevels));
bcrit_SNR = zeros(nDtype,length(nlevels));
Spower_MC = zeros(nDtype,length(nlevels),NITER);
eff_MC = zeros(nDtype,length(nlevels),NITER);
S2bias_MC = zeros(nDtype,NITER);
VBbias_MC = zeros(nDtype,NITER);
ratio=zeros(nDtype, NITER);
close all;

for nl=1:length(nlevels)
    fprintf('\rNoise level ...  %f ', nlevels(nl));
    delta = nlevels(nl);
    simp = zeros(NITER, tlen/2);
    run = zeros(NITER,tlen-1);
    sinterp=zeros(NITER,tlen);
    orig=zeros(NITER,tlen);
    sincc=zeros(NITER,tlen);

    VB_raw = zeros(NITER,2);
    VB_simp = zeros(NITER,2);
    VB_run = zeros(NITER,2);
    VB_sur = zeros(NITER,2);
    VB_sinc = zeros(NITER,2);

   
    %%% do the estimation for each case

    for iter=1:NITER
        fprintf('\rSimulation number: %d noise = %f ', iter, nlevels(nl));


        if (doNewMatrix==1)
            fprintf(' Making new design matrix ...')
            % Make the design matrix (main effect)
            X = zeros(tlen,1);
            %X(round([18:10:360]/1.4))=1;


            % 1=blocked, 2=ER Fixed, 3=ER randomized
            switch abs(DesignType)
                case 1
                    isi = [1:30 60:90 120:150 180:210 240:tlen];
                case 2
                    isi = [1:20:600 2:20:600]/1.4;
                case 3
                    isi = rand(30,1)*12 + 5;
                    isi = round(cumsum(isi));
                    isi = isi(find(isi<tlen));
            end
            X(round(isi))=1;
            X = conv(X,spm_hrf(TR));
            X = X(1:tlen);
            X = X/max(X);
            BOLD=X; %save a copy of BOLD so we can use X for other things
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
	      if noBOLD   %so, this creates a model without a BOLD regressor
	      % modulation of the effect
	      mm = ones(size(X));
	      mm(1:2:end)=-1;
	      X=X.*mm/2;
	      % baseline perfusion: regressor for controls and regressor for tag:
	      Xc = zeros(tlen,1);
	      Xt = zeros(tlen,1);
	      Xc(2:2:end)=1;
          Xc=Xc-.5;
	      Xt(1:2:end)=1;                  
	      Xb = ones(size(Xc));    % intercept
          X = [Xb Xc  X];
	    else
	      %Build BOLD design mtx
	      mm = ones(size(X));
	      mm(1:2:end)=-1;
	      X=X.*mm/2;
	      Xc = zeros(tlen,1);
	      Xt = zeros(tlen,1);
	      Xc(2:2:end)=1;
          Xc=Xc-0.5;
	      Xt(1:2:end)=1;
 	      Xb = ones(size(X));
          X = [Xb Xc BOLD X];
          end
          
        else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %    load /Users/hernan/matlab/flow/noise/ASLX.mat
            fprintf(' Reading design matrix ...')
	    if isASL
	      load ASLrER.mat
	    else
	      load BOLDrER.mat
	    end
        end
        

	if noBOLD
	  Bint = logical([0 0 1]);  % Beta's of interest
    else
	  Bint = logical([0 0 0 1]);    % Beta's of interest
    end


    Vo = WKfun2('AR+WN',rho,VarAR,VarWN,tlen);  %This is the new and improved AR+WN
    V=Vo;  %THis is a reduntancy We're no longer changing SNR by altering the variance, we're changing the signal below
    W = WKfun2('mkW',[],V);   % W = V^(-1/2) to make inversions easier down the line...
    t=0:tlen-1;


    % method 1: no differencing
    D=D1; 
    DX = D*X;
    V2=D*V*D';
    V_d=V2/V2(1,1);        %this is the correlation
    s2_d=V2(1,1);          %this is the variance
        
	Bi = Bint;
	zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = []; 
    np_D1=length(Bi);  %This is used later to compute DF
    
	gls_VB_raw = Bi*inv(X'*inv(V)*X)*Bi';  % recall that W = V^(-1/2) only used w D1
    
	ols_VB_raw = s2_d*Bi*pinv(DX)*V_d*(pinv(DX))'*Bi';  
    
    
	ols_S2bias_raw = (trace(V_d) - trace(DX'*V_d*DX*(DX'*DX)^(-1))) /(size(DX,1)-size(DX,2));
   
    ols_VBbias_raw = ols_S2bias_raw *( Bi*inv((DX)'*(DX))*Bi')/( Bi*pinv(DX)*V_d*(pinv(DX))'*Bi' ) ; 
    
    c1=[0 0 1 0]; c2=[0 0 1 0];
    rat_D1=s2_d*c1*pinv(DX)*V_d*(pinv(DX))'*c1'/(c2*pinv(X)*V*(pinv(X))'*c2');
    
    % 	ols_VBbias_raw = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ( Bi*ols_VB_raw*Bi' ) * ols_S2bias_raw ; 
    
    % method 2: simple subtraction

    D=D2;
    DX = D*X; 
    V2=D*V*D';       %note that V has the variance and correlation combined.
    V_d=V2/V2(1,1);  %correlation
    s2_d=V2(1,1);    %variance
	Bi = Bint;
	zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = []; 
    np_D2=length(Bi);
	
	ols_VB_simp = s2_d*Bi*pinv(DX)*V_d*(pinv(DX))'*Bi';
    
    %ols_VB_simp = s2_d*Bi*inv(DX'*DX)'*Bi';
    
    
    
	ols_S2bias_simp = (trace(V_d) - trace(DX'*V_d*DX*(DX'*DX)^(-1))) /(size(DX,1)-size(DX,2));
    
    ols_VBbias_simp =  ols_S2bias_simp * (Bi*inv(DX'*DX)*Bi')/( Bi*pinv(DX)*V_d*(pinv(DX))'*Bi' ); 
    
     c1=[0 1 0]; c2=[0 1 0 0];
    rat_D2=s2_d*c1*pinv(DX)*V_d*(pinv(DX))'*c1'/(c2*(pinv(X)*V*(pinv(X))')*c2');
%     ols_VBbias_simp = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ... 
% 	    ( Bi*ols_VB_simp*Bi' ) * ols_S2bias_simp; 

   
   % method 3: running subtraction

    D=D3;
    DX = D*X;
    V2=D*V*D';
    V_d=V2/V2(1,1);
    s2_d=V2(1,1);
	Bi = Bint;
	zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = []; 
    np_D3=length(Bi);
    
	ols_VB_run = s2_d*Bi*pinv(DX)*V_d*(pinv(DX))'*Bi';
    ols_S2bias_run = (trace(V_d) - trace(DX'*V_d*DX*(DX'*DX)^(-1))) /(size(DX,1)-size(DX,2));

    ols_VBbias_run =  ols_S2bias_run* ( Bi*inv(DX'*DX)*Bi')/(Bi*pinv(DX)*V_d*(pinv(DX))'*Bi' ) ;
    
    c1=[0 1 0]; c2=[0 1 0 0];
    rat_D3=s2_d*c1*pinv(DX)*V_d*(pinv(DX))'*c1'/(c2*(pinv(X)*V*(pinv(X))')*c2');
     
     % 	ols_VBbias_run = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ... 
     %  ( Bi*ols_VB_run*Bi' ) * ols_S2bias_run;  
     
     
     
     % method 4: surround subtraction

    D=D4;
    DX = D*X;
    V2=D*V*D';
    V_d=V2/V2(1,1);
    s2_d=V2(1,1);
    
	Bi = Bint;
	zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = []; 
    np_D4=length(Bi);
	
	ols_VB_sur = s2_d*Bi*pinv(DX)*V_d*(pinv(DX))'*Bi';
	ols_S2bias_sur = (trace(V_d) - trace(DX'*V_d*DX*(DX'*DX)^(-1)))/(size(DX,1)-size(DX,2));
    
    ols_VBbias_sur =ols_S2bias_sur*( Bi*inv(DX'*DX)*Bi')/ (Bi*pinv(DX)*V_d*(pinv(DX))'*Bi')  ; 
    
     c1=[0 1 0]; c2=[0 1 0 0];
    rat_D4=s2_d*c1*pinv(DX)*V_d*(pinv(DX))'*c1'/(c2*(pinv(X)*V*(pinv(X))')*c2');

    % 	ols_VBbias_sur = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ... 
% 	    ( Bi*ols_VB_sur*Bi' ) * ols_S2bias_sur; 

    
        
        
        % method 5: sinc subtraction

    D=D5;
    DX = D*X;
    V2=D*V*D';
    V_d=V2/mean(diag(V2));  %heteroscedastic, so I had to do it different mean(diag(V_d))=1
    s2_d=mean(diag(V2));
	Bi = Bint;
	zCol = find(all(abs(DX-0)<1e-6)); DX(:,zCol) = []; Bi(zCol) = []; 
    np_D5=length(Bi);
    
	
	ols_VB_sinc = s2_d*Bi*pinv(DX)*V_d*(pinv(DX))'*Bi';
    
	ols_S2bias_sinc = (trace(V_d) - trace(DX'*V_d*DX*(DX'*DX)^(-1)))/(size(DX,1)-size(DX,2));

    ols_VBbias_sinc = ols_S2bias_sinc*( Bi*inv(DX'*DX)*Bi') /( Bi*pinv(DX)*V_d*(pinv(DX))'*Bi' ) ; 
    
     c1=[0 0 1]; c2=[0 0 0 1];
    rat_D5=s2_d*c1*pinv(DX)*V_d*(pinv(DX))'*c1'/(c2*(pinv(X)*V*(pinv(X))')*c2');
    % 	ols_VBbias_sinc = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ... 
    % ( Bi*ols_VB_sinc*Bi' ) * ols_S2bias_sinc; 

    nparams=...
        [np_D1;
         np_D2;
         np_D3;
         np_D4;
         np_D5;
         np_D1;
         ];
     
	var_Bhat = ...
	    [ols_VB_raw;
	     ols_VB_simp;
	     ols_VB_run;
	     ols_VB_sur;
	     ols_VB_sinc;
	     gls_VB_raw;
	    ];
	
	S2bias = ...
	    [ols_S2bias_raw;
	     ols_S2bias_simp;
	     ols_S2bias_run;
	     ols_S2bias_sur;
	     ols_S2bias_sinc;
	     1
	    ];
	
	VBbias = ...
	    [ols_VBbias_raw;
	     ols_VBbias_simp;
	     ols_VBbias_run;
	     ols_VBbias_sur;
	     ols_VBbias_sinc;
	     1
	    ];
    
     ratio=...
         [rat_D1;
          rat_D2;
          rat_D3;
          rat_D4;
          rat_D5;
          1];

        eff = 1./var_Bhat;
        eff_rel = eff ./ eff(end);


        % Power calculation given mean, std.dev, and known effect size
        % we'll test the power of estimation of the second B param. in all
        % three types of differencing
       
        q = repmat(NaN,nDtype,1);
        alpha = 0.0001;
              
        
        for Dtype=1:nDtype    % repeat for each differencing method
            
            df=length(X)-length(nparams(Dtype));
            
            effect_size=delta/sqrt(var_Bhat(Dtype));
            
            sigma2_B = var_Bhat(Dtype);

            % find which value of beta_hat corresponds to the
            % significance level alpha
            % in the null dustribution:
            tcrit =  spm_invTcdf(1-alpha, df);
            bcrit = tcrit * sqrt(sigma2_B);
            q(Dtype) = spm_Ncdf(bcrit, effect_size, sigma2_B);

            bcrit_SNR(Dtype,nl) = bcrit;
        end

        Spower = 1-q;

        % and of course over several noise levels:
       
        Spower_SNR (:,nl) = Spower;
        eff_SNR (:,nl) = eff;

        % in the case of a monte carlo simulation, we'll average the power
        % and efficiency over iterations.  This means that we are not
        % keeping track of the noise level any more!
        Spower_MC (:,nl,iter) = Spower;
        eff_MC (:,nl,iter) = eff;
	S2bias_MC (:,nl,iter) = S2bias;
	VBbias_MC (:,nl,iter) = VBbias;

        % make a little table
        %     fprintf('\n Power  Efficiency   \n');
        %     table = [Spower(:,3)  eff(:,3)  ]
        %     table ./ repmat(table(1,:),5,1)

    end


    table = [Spower eff Spower./Spower(end) eff./eff(end) S2bias VBbias];
    fprintf('\n    Power      Eff      RelPow    RelEff    S2bias    VBbias\n');
    disp(table)
end
