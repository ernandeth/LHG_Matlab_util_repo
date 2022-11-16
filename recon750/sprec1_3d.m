function sprec1_3d_par(strFile,varargin)

% Deal with arguments
[args,scaninfo,kinfo] = rec_setup1(strFile,varargin{:});

% kludge: read in the additional (fake) frame. 
oddframes = 0;
if rem(scaninfo.nphases,2)>0
    scaninfo.nphases = scaninfo.nphases+1
    oddframes  = 1;
end

if args.iph == 0
    iPhase = 1:scaninfo.nphases;
else
    iPhase = args.iph;
end

parfor t = 1:scaninfo.nphases
%for t = iPhase
    
    % Change args for each TR 
    argsNew = args;
    argsNew.iph = t;  

    str_file = do_sprec1_3d(argsNew,scaninfo,kinfo);
    

end

% delete(gcp);

% Combine into one 4D .nii

vfiles = dir('vol_e*.nii');

% kludge: remove the last image if odd frames
if oddframes
    vfiles = vfiles(1:end-1);
end

if ~isempty(vfiles)
    strDate = datestr(now,'yyyy-mm-dd_HHMMSS');
    read_nii_series('vol_e*',[strDate '.nii'] );
    !rm -f vol_e*nii
    str= ['!mv ' strDate '.nii ' vfiles(end).name]
    eval(str)
end

pvfiles = dir('p_vol_e*.nii');
if oddframes
    pvfiles = pvfiles(1:end-1);
end

if ~isempty(pvfiles) && args.complex
    strDate = datestr(now,'yyyy-mm-dd_HHMMSS');
    read_nii_series('p_vol_e*',[strDate '.nii'] );
    !rm -f p_vol_e*nii
    str= ['!mv ' strDate '.nii ' pvfiles(end).name]
    eval(str)
end
