function [resp_phases,card_phases,resp,card,respnew,tr_samp]=get_phys_phase(physfile,tr,nslices,npr,nextra,nframes,samprate);
%
% Usage: [resp_phases,card_phases,resp,card,respnew,tr_samp]=get_phys_phase(physfile,tr,nslices,npr,nextra,nframes,samprate);
%
%  This extracts the physiological phases from the recorded physiological data, which can then be used with the function phys_correct
%  to implement physiological correction for fMRI data.
%
%  Inputs       physfile                  -  the externally recorded physio data; it should be in column text format, with no text headers.
%               tr                        -  the repetition time of the scan, in seconds.
%               nslices                   -  number of scans.
%               npr                       -  number of shots in the scan.
%               nextra                    -  number of disdaqs collected in front of the scan.
%               nframes                   -  number of temporal frames in the scan.
%               samprate                  -  A/D sampling rate used, in Hz (should be set at 100 Hz).
%
%  Outputs      resp_phases, card_phases  -  ordered relative physiological phases (between 0 and 1), with size [nslices npr nframes]
%               resp, card                -  recorded physio data from physfile, use to check if input channels are correct. 
%               respnew                   -  filtered respiratory data, use to check if filter range is correct.
%               tr_samp                   -  locations of each tr acquisition.
%
%%%%%%%%%%%%%%%%%%%%%

%% Make sure ordering is correct

time=physfile(:,1);         %  Time vector, not used, but should be monotonically increasing
resp=physfile(:,2);         %  Respiration waveform
card=physfile(:,3);         %  Cardiac waveform (only peak locations from GE)

xgrd=physfile(:,4);         %  Gradient waveform (used for timing info, x or y is fine for spiral, may matter for other sequences)

%%%%%%%%%%%%%%%%%%%%%

%%% Smoothing of respiratory waveform to frequency band of interest


sampfreq=1/samprate;

low_cutoff=floor(0.1*sampfreq*length(resp));
high_cutoff=ceil(0.5*sampfreq*length(resp));

bnd_pss_width=high_cutoff-low_cutoff;

filt_resp=[zeros(1,low_cutoff) ones(1,bnd_pss_width) zeros(1,(length(resp)-(2*(low_cutoff+bnd_pss_width)))) ones(1,bnd_pss_width) zeros(1,low_cutoff)];

sg=fft(resp.').*filt_resp;

respnew=real(ifft(sg))+mean(resp);


%%% Extract respiratory phase

x=1:length(resp);

s=[x.' respnew.'];

big_resp_phs=res_phs_new(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Extract cardiac phase

ctst=(card<0);

s=[x.' ctst];

big_card_phs=res_phs_new2(s);

%%%%%%%%%%%%%%%%%%%%%%%%%5

%%% Extract tr locations

lst=find(abs(xgrd)>.1);     %  Threshold for gradient waveform being on, might want to verify

x_strt=lst(1);
x_end=lst(length(lst));

skip=tr/(sampfreq);

nviews=nslices;

%tr_samp=x_strt:skip/nviews:length(xgrd);
tr_samp=x_strt:skip/nviews:x_end;
tr_samp=round(tr_samp(1:(nviews*nextra + nslices*npr*nframes)));
keyboard;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Sample respiratory and cardiac at each TR

resp_samp=big_resp_phs(tr_samp);

card_samp=big_card_phs(tr_samp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Remove disdaqs and reshape

resp_phases=reshape(resp_samp(((nviews*nextra)+1):(length(resp_samp))),[nslices npr nframes]);

card_phases=reshape(card_samp(((nviews*nextra)+1):(length(resp_samp))),[nslices npr nframes]);
