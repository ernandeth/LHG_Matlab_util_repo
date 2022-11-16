function [Yres,df]=Blind_phys_corr_PCA(Y,TR,istime);
% -----------------------------------
% Usage: [Yres,df]=Blind_phys_corr_PCA(Y);
% -----------------------------------
% Does supervised physio correction.
% -----------------------------------
% Input:  Y:      TxM fMRI Data matrix, T is 
%                 #time points, M is #voxels
%         TR:     Sampling time in sec
%         istime: if == 1 then time domain o.w. freq. domain
% -----------------------------------
% Output:  Yres:    fMRI residual TxM data matrix
%          df:      degrees of freedom 
% ----------------------------------
% Magnus Orn Ulfarsson, 2007.
% -----------------------------------
if(nargin==2)
    istime=1;
end
[T,M]=size(Y);
mu=mean(Y)';
Y=Y-ones(T,1)*mu';
df=T-1;
t=linspace(0,T*TR,T)';
[r]=laplace_pca(Y);
[G,s2,l]=nPCA(Y,0,r);
if(r>64), r=64; end
P=Y*G;
if(istime==0)
    for j=1:r,
        fs=1/TR;
        [Fy,f]=periodogram_mu(P(:,j),fs,[0,fs/2]);
        Pfreq(:,j)=Fy;
    end
end


numFigs=ceil(r/16);
numPlots_lf=mod(r,16);
num=1;
for i=1:numFigs,
    figure(i);
    if(i==numFigs), 
        for j=1:numPlots_lf, 
            subplot(4,4,j)
            if(istime==1), 
                plot(t,P(:,j)); xlabel('t [sec]'); title(num2str(num)); num=num+1;
            else
                plot(f,Pfreq(:,j)); xlabel('f [Hz]'); title(num2str(num)); num=num+1;
            end           
        end
    else
         for j=1:16, 
            subplot(4,4,j) 
            if(istime==1), 
                plot(t,P(:,j)); xlabel('t [sec]'); title(num2str(num)); num=num+1;
            else
                plot(f,Pfreq(:,j)); xlabel('f [Hz]'); title(num2str(num)); num=num+1;
            end           
         end
    end
end
disp(['Type in the matrix Ind which PCs correspond to Physio noise'])
disp(['For example: Ind=[1 2 4];']);
disp(['Then type return']);
keyboard
reff=length(Ind);
df=df-reff;
R=eye(r)+s2*inv(diag(l));
Reff=eye(reff)+s2*inv(diag(l(Ind)));
Yres = Y*G*inv(R)*G'-Y*G(:,Ind)*inv(Reff)*G(:,Ind)';
Yres=Yres+ones(T,1)*mu';




