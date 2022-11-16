function [DX,sf] = tor_make_deconv_mtx2(sf,tp,eres)
% function [DX,sf] = tor_make_deconv_mtx(sf,tp,eres)
%   sf: cell array of stick functions, one per condition
%       all sf cells should be of the same length
%
%   tp: number of timepoints to estimate in hrf deconvolution matrix
%   eres: timebins in sf array for each TR
%
%   DX: deconvolution matrix
%       estimates O.tp time points for each condition
%       Time resolution is in TRs
%   
%   sf: stick function resampled at TR
%
%   No parametric modulation of sf's allowed.
%
% Tor Wager, 10/20/01

shiftElements = eres;
% each time point is timeRes TRs.

% -------------------------------------------------------------------
% * downsample sf to number of TRs
% -------------------------------------------------------------------
numtrs = length(sf{1}) ./ eres;
myzeros = zeros(numtrs,1);

for i = 1:length(sf)
    Snumtrs = length(sf{i}) ./ eres;
    if Snumtrs ~= round(Snumtrs), warning(['sf{ ' num2str(i) '}: length not evenly divisible by eres.']),end
    if numtrs ~= Snumtrs, warning(['sf{ ' num2str(i) '}: different length than sf{1}.']),end

    inums = find(sf{i} > 0);
    inums = inums ./ eres;      % convert to TRs
    inums = ceil(inums);       	% nearest TR
    inums(inums == 0) = 1;		% first possible TR is 1st element
    sf{i} = myzeros;
    sf{i}(inums) = 1;           % always use 1 for sf
end

% -------------------------------------------------------------------
% * make deconvolution matrix DX
% -------------------------------------------------------------------

index = 1;
for i = 1:size(sf,2)
    DX(:,index) = sf{i};
    index = index + 1;
    inums = find(sf{i} == 1);
    
    for j = 2:tp
        inums = inums + 1;
        reg = myzeros;
        reg(inums) = 1;
        reg = reg(1:numtrs);
        DX(:,index) = reg;
        index  = index + 1;
    end
    
end

% -------------------------------------------------------------------
% * add intercept
% -------------------------------------------------------------------
DX(:,end+1) = 1;


return