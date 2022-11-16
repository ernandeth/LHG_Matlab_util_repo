function dat = getdespikeddat(phnum,slnum,scaninfo,strRawFile)

strCurDir = pwd;
[strPath,strName] = fileparts(strRawFile);

if ~isempty(strPath)
    cd(strPath);
end

strPat = sprintf('f*%s.data',strName);
% s = dir('f*crap.data');
s = dir(strPat);
casFiles = {s(:).name}';
casFiles = sort(casFiles);

ndat = scaninfo.ndat;
% nphases = scaninfo.nphases; % check for even/odd thing
nphases = scaninfo.nphmult;
ncoils = scaninfo.ncoils;

strFile = casFiles{slnum+1};
fidSl = fopen(strFile,'rb');

% offs = phnum * 4 * ndat * ncoils;
% fseek(fidSl,offs,'bof');

dat = zeros(ndat,ncoils);
for coilnum = 1:ncoils
    coilOff = (coilnum-1) * 4 * ndat * nphases;
    phaseOff = phnum * 4 * ndat;
    fseek(fidSl,coilOff+phaseOff,'bof');
    ktmp = fread(fidSl,2*ndat,'int16');
    dat(:,coilnum) = ktmp(1:2:end) + i.*ktmp(2:2:end);
end
fclose(fidSl);

cd(strCurDir);