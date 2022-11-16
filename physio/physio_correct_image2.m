function [fid_new,sldat,sldat_new]=physio_correct_image2(data,order,resp_flag,card_flag,resp_phases,card_phases,strt_time);
%
% Usage:  [fid_new,sldat,sldat_new]=physio_correct(fname,order,resp_flag,card_flag,resp_phases,card_phases,strt_time) 
%
%  This implements physiological noise correction for k-space fMRI data, along with detrending and outlier removal.  
%
%  Inputs       data                      -  original image space data, with size [xres yres nphases nslices]
%               order                     -  order of polynomial to use for detrending.
%               resp_flag, card_flag      -  if 1, do physio correction, if 0, do not
%               resp_phases, card_phases  -  ordered relative physiological phases (between 0 and 1), with size [nslices 1 nframes]
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

RAWHEADERSIZE=39984;   % Header size for GE Signa LX, usually changes with software updates


% Defining files

xres=size(data,1);
yres=size(data,2);
nphases=size(data,3);
nslices=size(data,4);

data_nw=zeros(size(data));


for k=1:nslices,  %nslices

		disp(['Initiating Analysis...']);

        
        %% Form regressor matrix for respiration
        r_phases=squeeze(resp_phases(k,1,strt_time:size(resp_phases,3)));
        
		A=[ones([length(r_phases) 1]) cos(2*1*pi*r_phases) sin(2*1*pi*r_phases) cos(2*2*pi*r_phases) sin(2*2*pi*r_phases)];
		D=ones(size(A)); M=mean(A); M(1)=0;
		for m=1:size(A,2), D(:,m)=M(m)*D(:,m); end;
		A=A-D; %size(A)
        
        
        %% Form regressor matrix for cardiac
        c_phases=squeeze(card_phases(k,1,strt_time:size(card_phases,3)));
        
		B=[ones([length(c_phases) 1]) cos(2*1*pi*c_phases) sin(2*1*pi*c_phases) cos(2*2*pi*c_phases) sin(2*2*pi*c_phases)];
		F=ones(size(B)); N=mean(B); N(1)=0;
		for m=1:size(B,2), F(:,m)=N(m)*F(:,m); end;
		B=B-F; %size(A)
        
        
  		for x=1:xres
            for y=1:yres
			
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%	Ignore initial timepoints, if needed
 
            tcorg=squeeze(data(x,y,:,k)).';
            
            tco=tcorg(strt_time:length(tcorg));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Detrend  (mandatory)
            
    	    tr_val=polyfit(1:length(tco),tco,order);
    		
            tr_fit=polyval(tr_val,1:length(tco));
            
    		det_tc=tco-tr_fit+mean(tco);
    		
            tcnew(strt_time:length(tcorg))=det_tc;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%  Correct for respiration
            
            if (resp_flag)
                
                
                coeff_resp=A\det_tc.'; %size(coeff_re)
                
                resp_tc=det_tc - (A*coeff_resp).' + mean(det_tc);
                
                tcnew(strt_time:length(tcorg))=resp_tc;
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
  
            data_nw(x,y,:,k)=tcnew;    
          
        end;  %xres
    end; %yres

disp(['Slice = ',int2str(k)]);    
    
end;  % slice

disp('Physio correction done!!')