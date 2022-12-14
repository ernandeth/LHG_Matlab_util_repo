% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Robert C. Welsh
% Copyright 2005
%
%
% A routine to bandpass filter data
%
%     theData     = theData(space,time) (this is the 
%                   standard format being used in this SOM
%                   implementation).
%
%     sample      = sample period (TR in fmri language)
%
%     lowFreq     = low frequency cutoff
%
%     highFreq    = high frequency cutoff
%
%     gentle      = 0 - no just a hard cut, 1 - yes.
%                   gentle is a curve of 0.1, .5, .9, 1.0 on the rising and
%                   falling edges except for DC.
%                 = 0 to use fir2 when using the Signal Processing ToolBox
%                   1 to use firpm 
%
%     padding      = padding for filtering. (default is 10)
%
%     whichFilter = 0 try to use TBX, else use SOM_Filter_FFT
%
% function [results, b, filtParms] = SOM_Filter(theData,sample,lowFreq,highFreq,gentle,padding,whichFilter)
%
%
% If the signal toolbox is present then use that to do the
% filtering.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [results, b, filtParms] = SOM_Filter(theData,sample,lowFreq,highFreq,gentle,padding,whichFilter)

% Save input for posterity.

filtParms.sample      = sample;
filtParms.lowFreq     = lowFreq;
filtParms.highFreq    = highFreq;
filtParms.gentle      = gentle;
filtParms.padding     = padding;
filtParms.whichFilter = whichFilter;

% make sure padding exists.

if exist('padding') == 0
  fprintf('Default padding being set to 10.\n');
  padding = 10;
end

% make sure "whichFilter" is existing.

if exist('whichFilter') == 0
  fprintf('Default to trying to use Signal Processing Toolbox.\n');
  whichFilter = 0;
end

% If padding is on then we need to pad the beginning and the end 
% by "padding" number of average images. These are then 
% removed later.

% Make sure padding is integer and > 0.

padding = floor(padding);

if padding < 0
  padding = 0;
end

theDataMean = mean(theData,2);

theDataPad  = repmat(theDataMean,1,padding);

theData = [theDataPad theData theDataPad];

% Determine if firpm is present and filtfilt
% Get the fft of the data.

%if whichFilter ~= 0
  fprintf('Using handmade fft filter method.\n');
  [results b] = SOM_Filter_FFT(theData,sample,lowFreq,highFreq,gentle);
%else
%  if exist('firpm.m') == 2 & exist('filtfilt.m') == 2 & whichFilter == 0
%  fprintf('Using Signal Processing ToolBox\n');
%  [results b] = SOM_Filter_SIGTBX(theData,sample,lowFreq,highFreq,gentle);
%end

% Now trim out the padding;

results = results(:,1+padding:end-padding);

return

%
% Using the Signal Processing ToolBox
%

function [results, b] = SOM_Filter_SIGTBX(theData,sample,lowFreq,highFreq,firOpt);

%Determine the Nyquist criterion

nyquist = 1/sample/2;

% How many time points.

N = size(theData,2);

nFilt = floor(N/3)-1;

if nFilt < 1
  fprintf('Error - insufficient data to actually filter.\n');
  results = [];
  b = [];
  return
end

% Need to put in protection for 
%     1) lowFreq <=0  and/or
%     2) highFreq >=nyquist.

if lowFreq > 0 
    fir2Freq = [0 lowFreq/nyquist*.999 lowFreq/nyquist ];
    fir2coef = [0              0              1        ];
else
    fir2Freq = [0 ];
    fir2coef = [1 ];
end

if highFreq < nyquist 
    fir2Freq = [fir2Freq highFreq/nyquist highFreq/nyquist*1.001 1.0 ];
    fir2coef = [fir2coef         1                 0              0  ];
else
    fir2Freq = [fir2Freq 1.0 ];
    fir2coef = [fir2coef   1 ];
end

if firOpt == 0
  b = fir2(nFilt,[0 lowFreq/nyquist*.999 lowFreq/nyquist highFreq/nyquist ...
		  highFreq/nyquist*1.001 1.0],[0 0 1 1 0 0]);
else
  b = firpm(nFilt,[0 lowFreq/nyquist*.999 lowFreq/nyquist highFreq/nyquist ...
		  highFreq/nyquist*1.001 1.0],[0 0 1 1 0 0]);
end

% Filter the data now, forwards and backwards, using 
% the transpose of data as time flows down, not across.

results = (filtfilt(b,1,theData'))';

return


%
% Using the FFT Method. Problem with Phase?
%

function [results, b] = SOM_Filter_FFT(theData,sample,lowFreq,highFreq,gentle)


%Determine the Nyquist criterion

nyquist = 1/sample/2;

% How big is our sample.

N = size(theData,2);

% Make a frequency baseline.

deltaF = nyquist/(floor(N/2)-1);

freq   = (-floor(N/2):floor(N/2)-1)*deltaF;

% Get the fft of the data.

ffttheData = fftshift(fft(theData,[],2),2);

% Find what to remove.

fftMask = zeros(size(theData));
maskH  = zeros(size(theData));
maskL  = zeros(size(theData));

% Need to put in protection for 
%     1) lowFreq <=0  and/or
%     2) highFreq >=nyquist.

ifrqHi = find(abs(freq)<=highFreq);
ifrqLo = find(abs(freq)>=lowFreq);

maskH(:,ifrqHi) = 1;
maskL(:,ifrqLo) = 1;

fftMask = maskH.*maskL;

% Find the transition points and roll a little.

iup = find(diff(fftMask(1,:))>0);
idn = find(diff(fftMask(1,:))<0);

if gentle ~= 0
    fftMask(:,iup-1) = .1;
    fftMask(:,iup) = .5;
    fftMask(:,iup+1) = .9;
    fftMask(:,idn+2) = .1;
    fftMask(:,idn+1) = .5;
    fftMask(:,idn) = .9;
end

% DC Component is saved.

dci = floor(N/2)+1;

fftMask(:,dci) = 1;

% Filter and return the results - returning the real component.

results = real(ifft(ifftshift(fftMask.*ffttheData,2),[],2));

b = fftMask;

return

%
% All done.
%

