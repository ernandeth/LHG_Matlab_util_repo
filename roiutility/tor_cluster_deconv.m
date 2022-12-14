function [dx,xY,b,hrf] = tor_cluster_deconv(varargin)
% function [dx,xY,b,hrf] = tor_cluster_deconv(varargin)
%
% arguments: string followed by value
% 'data', 		y vector to fit DX to
% 'tp', 		time points (how many TRs) to deconvolve
% 'samprate', 	 	sampling rate of stick function in Sess
% 'columns', 	 	which columns of stick function are of interest in deconvolution
% 'Sess', 		Sess structure from SPM.mat
% 'DX', 		deconvolution matrix, if it exists.

% -------------------------------------------------------------------
% * set up input arguments
% -------------------------------------------------------------------

dx = [];
y = [];
tp = [];
eres = [];
colsOfInterest = [];
Sess = [];
DX = [];
doplot = 0;

for i = 1:length(varargin)
	if isstr(varargin{i})
		switch varargin{i}
		case 'data', 		y = varargin{i+1};
		case 'tp', 			tp = varargin{i+1};
		case 'samprate', 	eres = varargin{i+1};
		case 'columns', 	colsOfInterest = varargin{i+1};
		case 'Sess', 		Sess = varargin{i+1};
		case 'DX', 			DX = varargin{i+1};
		case 'plot',		doplot = 1;
		end % end switch
	end
end


% -------------------------------------------------------------------
% * make deconvolution matrix if necessary
% -------------------------------------------------------------------

if isempty(DX)
	if isempty(Sess) | isempty(tp),	error('Must enter either deconvolution matrix DX or Sess struct and tp'), end
	if isempty(eres), warning('No sampling rate for stick function entered: using SPM default of 16 timepoints per TR'),end

	[DX,sf,hires_sf] = makeDX(Sess,tp,eres,colsOfInterest);
	
    dx.hires_sf = hires_sf;
	dx.DX = DX;
	dx.sf = sf;
end


% -------------------------------------------------------------------
% * fit to y if there's data
% -------------------------------------------------------------------

if ~isempty(y)
	if length(y) ~= size(DX,1),
		whos y, whos DX, error('DX and y are different lengths.')
	end

	b = pinv(DX) * y;
	xY.b = b;

	% -------------------------------------------------------------------
	% * break betas into hrf estimates if possible
	% -------------------------------------------------------------------
	if ~isempty(tp)
		index = 1;
		for i = 1:tp:length(b)-1
			hrf{index} = b(i:i+tp-1);
			index = index + 1;
		end
		xY.hrf = hrf;
	end
	

end


% -------------------------------------------------------------------
% * plot, if all required variables are there
% -------------------------------------------------------------------
if doplot == 1 & exist('hrf') == 1

	mycolors = {'ro-' 'go-' 'bo-' 'ko-' 'mo-' 'r^-' 'g^-' 'b^-' 'k^-' 'm^-' 'y^-'};
	cla, hold on, grid on
	title('Hemodynamic responses')
	whos hrf
	for i = 1:length(hrf)
		plot(hrf{i},mycolors{i},'LineWidth',2)
		myleg{i} = ['Condition ' num2str(i)];
	end
	legend(myleg)
end


return