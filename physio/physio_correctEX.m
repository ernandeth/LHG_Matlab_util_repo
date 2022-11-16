function [fid_new,sldat,sldat_new]=physio_correctEX(fname,order,resp_flag,card_flag,resp_phases,card_phases,strt_time, Nframes);
%
% Usage:  [fid_new,sldat,sldat_new]=physio_correct(fname,order,resp_flag,card_flag,resp_phases,card_phases,strt_time) 
%
%  This implements physiological noise correction for k-space fMRI data, along with detrending and outlier removal.  
%
%  Inputs       fname                     -  name of raw Pfile to correct (e.g. 'P10001.7')
%               order                     -  order of polynomial to use for detrending.
%               resp_flag, card_flag      -  if 1, do physio correction, if 0, do not
%               resp_phases, card_phases  -  ordered relative physiological phases (between 0 and 1), with size [nslices npr nframes]
%                                             these can be generated with function get_phys_phase
%               strt_time                 -  frame index to start correction, the range strt_time:end will be corrected, with the range
%                                             1:(strt_time-1) being untouched.
%  Outputs      
%               fid_new                   -  file_identifier for the corrected Pfile, if it is -1, there is an error in opening the file
%               sldat                     -  the original slice data (when returned, will be the last slice)
%               sldat_new                 -  the corrected slice data (when returned, will be the last slice)
%%%%%%%%%%%%%%%%%%%%%

disp(['Det. Order:',int2str(order),' Resp.: ',int2str(resp_flag),' Card.: ',int2str(card_flag)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RAWHEADERSIZE=39984;   % Header size for GE Signa LX, usually changes with software updates
RAWHEADERSIZE=60464;   % Header size for GE Signa LX, usually changes with software updates
RAWHEADERSIZE=61464;   % Header size for GE Excite, usually changes with software updates

% Reading header of original pfile, spiral format
%[info1,info2,info3]=fidread2(fname);
[info1,info2,info3]=fidread2EX(fname);

ndat=info1(9)
nslices=info1(3)
npr=info3(5)
%nph=info3(2)
nph=Nframes
nphmult=info3(11)
nphases=nph*nphmult
frsize=2*npr*ndat

odd_flag=mod(nph*npr,2);    %should be right, check

disp(['odd_flag=',int2str(odd_flag)]);

% Opening files

if (resp_flag)
    if(card_flag)
        fname_new=[fname,'.d3rc'];
    else
        fname_new=[fname,'.d3r'];
    end;
else
    fname_new=[fname,'.d3'];
end;


fid_old=fopen(fname,'r','ieee-le');
fid_new=fopen(fname_new,'w','ieee-le');

if fid_old<3, error('Could not open old file...'); end;
if fid_new<3, error('Could not open new file...'); end;

% Write header to corrected pfile
[header_buff,ver]=fread(fid_old,RAWHEADERSIZE,'uint8');
ver=fwrite(fid_new,header_buff,'uint8');
if ver~=RAWHEADERSIZE, error('Could not write header!'); end;

for k=1:nslices,  %nslices
     
     %% Read and write baseline data
     bsln=fread(fid_old,[2*ndat],'short');
     fwrite(fid_new,bsln,'short');
     
     
     %% Read in slice data
     [sldat,slcnt]=fread(fid_old,[2*npr*ndat*nphases],'short');
     
     if (slcnt~=2*npr*ndat*nphases)
         error(['Incorrect amount of data! ',int2str(slcnt),'/',int2str(2*npr*ndat*nphases)]);
     end;
     
     sldat_nw=zeros(size(sldat));
     
	for n=1:npr,

		disp(['Initiating Analysis...']);

        
        %% Form regressor matrix for respiration
        r_phases=squeeze(resp_phases(k,n,strt_time:size(resp_phases,3)));
        
		A=[ones([length(r_phases) 1]) cos(2*1*pi*r_phases) sin(2*1*pi*r_phases) cos(2*2*pi*r_phases) sin(2*2*pi*r_phases)];
		D=ones(size(A)); M=mean(A); M(1)=0;
		for m=1:size(A,2), D(:,m)=M(m)*D(:,m); end;
		A=A-D; %size(A)
        
        
        %% Form regressor matrix for cardiac
        c_phases=squeeze(card_phases(k,n,strt_time:size(card_phases,3)));
        
		B=[ones([length(c_phases) 1]) cos(2*1*pi*c_phases) sin(2*1*pi*c_phases) cos(2*2*pi*c_phases) sin(2*2*pi*c_phases)];
		F=ones(size(B)); N=mean(B); N(1)=0;
		for m=1:size(B,2), F(:,m)=N(m)*F(:,m); end;
		B=B-F; %size(A)
        
        
  		for kpt=1:2*ndat,
			
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%	Ignore initial timepoints, if needed

            tcorg = (sldat((kpt+(n-1)*2*ndat):frsize:length(sldat))).';
            
            tco = tcorg(strt_time:length(tcorg));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Detrend  (mandatory)
            
    	    tr_val = polyfit(1:length(tco),tco,order);
    		
            tr_fit = polyval(tr_val,1:length(tco));
            
            det_tc = tco-tr_fit+mean(tco);
    		
            tcnew(strt_time:length(tcorg)) = det_tc;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%  Correct for respiration
            
            if (resp_flag)
                
                coeff_resp = A\det_tc.'; %size(coeff_re)
                
                resp_tc = det_tc - (A*coeff_resp).' + mean(det_tc);
                
                tcnew(strt_time:length(tcorg)) = resp_tc;
            end; % resp_flag
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %  Correct for cardiac
            
            if (card_flag)
                
                coeff_card=B\resp_tc.'; %size(coeff_re)
                
                card_tc=resp_tc - (B*coeff_card).' + mean(resp_tc);

                tcnew(strt_time:length(tcorg))=card_tc;
                
            end; %card_flag
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Remove outliers
            
%			mask_out=(abs(resp_tc-mean(resp_tc)) > 2.5*std(resp_tc));

 % 			outlr_tc=(1-mask_out).*resp_tc  + mask_out*mean(resp_tc);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Compose new data
            
            tcnew(1:(strt_time-1))=tcorg(1:(strt_time-1));
        
            sldat_nw((kpt+(n-1)*2*ndat):frsize:length(sldat))=tcnew;    
            
            
        end;  %ndat
        
        disp(['Slice = ',int2str(k),'Pr=',int2str(n)]);
        
    end; %npr 

    %% Write corrected slice data
    
    fwrite(fid_new,sldat_nw,'short');
    
    if (odd_flag)
        bsln=fread(fid_old,[2*ndat],'short');
        fwrite(fid_new,bsln,'short');
    end;
    
end;  % slice

%% Close files

fclose(fid_old);

fclose(fid_new);

disp('Physio correction done!!')
