function [O,sub_clusters] = batch_multisubject_cluster_ui(varargin)
% function [O,sub_clusters] = batch_multisubject_cluster_ui(O)
% 
% O is option structure with fields:	
%	an_files				cell array, one cell per subject, single subject SPM.mat file with Xx
%	im_files				cell array, one cell per subject, padded matrix of all image file names
%	evt_onsets			cell array, one cell per subject, containing matrices of event onsets
%							one condition per column
%	snums					vector with numbers for subjects - saves as sub*_clusters
%	nsess					number of sessions in design
%
%	
%
% Tor Wager, 10/31/01

try		% try to do the whole thing. 

% ----------------------------------------------------------------------------------
% Process function inputs
% ----------------------------------------------------------------------------------
DX = [];
if nargin > 0,
	O = varargin{1};
	N = fieldnames(O);
else
	N = [];
end

% ----------------------------------------------------------------------------------
% User-entered inputs
% ----------------------------------------------------------------------------------

if ~any(strcmp(N,'snums')), O.snums = input('Enter vector of subject numbers to analyze: ');,end

if ~any(strcmp(N,'nsess')),O.nsess = input('How many sessions? ');,end

if ~any(strcmp(N,'nscans')), O.nscans = input('Enter number of scans per session: ');,end

if any(strcmp(N,'SPMdes'))
	SPMdes = O.SPMdes;
else
	SPMdes = input('Get clusters from SPM analysis or other mat file with variable ''clusters'' (1 or 0)?');
end
O.SPMdes = SPMdes;

if any(strcmp(N,'desLoad'))
	desLoad = O.desLoad;
else
	desLoad = 0;
	disp('Analyze using loaded design matrix or vectors of onset times?')
	disp('You still need to enter a design matrix for event times, to get other params.')
	desChoice = input('(Type design or events) ','s');
	switch desChoice
		case 'design', desLoad = 1;
		case 'events', desLoad = 2;
				if ~any(strcmp(N,'evt_onsets')),error('Must define evt_onsets in O input structure.'),end
		otherwise disp('Invalid choice.')
	end
end
O.desLoad = desLoad;

if ~any(strcmp(N,'HP')),O.HP = input('Enter HP filter length for analysis, in s, or blank for none: ');,end

if any(strcmp(N,'dodeconv'))
	dodeconv = O.dodeconv;
else
	dodeconv = input('Do deconvolution? 1 or 0: ');
end
O.dodeconv = dodeconv;

N = fieldnames(O);

% ----------------------------------------------------------------------------------
% Define or load clusters
% ----------------------------------------------------------------------------------
go = 0;

while ~go
	if SPMdes
		disp('Select SPM.mat file from which to derive clusters.')
		[SPM,VOL,xX,xCon] = spm_getSPM;
		[clusters] = tor_extract_rois([],SPM,VOL);
		go = 1;
	else
		if ~any(strcmp(N,'clDefMat'))
			O.clDefMat = spm_get([1],'*.mat','Select mat file containing variable ''clusters''',pwd,0);
		end
		try
			eval(['load ' O.clDefMat ' ''clusters'';'])	
			if isstruct(clusters),go = 1;,else disp('Variable clusters is not a structure.'),end
		catch
			disp('Error, or mat file does not contain variable named clusters.')
		end
	end
end

% ----------------------------------------------------------------------------------
% Crunch through subjects
% ----------------------------------------------------------------------------------


for mySubIndex = 1:length(O.snums)

	mySubject = num2str(O.snums(mySubIndex));
	disp(['Subject ' mySubject])
	disp('_____________________________________________________')

	% ----------------------------------------------------------------------------------
	% Get File Names
	% ----------------------------------------------------------------------------------

	if ~any(strcmp(N,'im_files')),O.im_files = [];,disp('No file names entered yet.'),end

	if length(O.im_files) < O.snums(mySubIndex)
		disp(['Subject ' num2str(O.snums(mySubIndex)) ' not entered yet.'])
		O.im_files{O.snums(mySubIndex)} = [];
	end
	
	if size(O.im_files{O.snums(mySubIndex)},1) < O.nscans .* O.nsess	
		disp(['Looking for ' num2str(O.nscans .* O.nsess) ' images, found ' num2str(size(O.im_files{O.snums(mySubIndex)},1))])
		O.im_files{O.snums(mySubIndex)} = spm_get([O.nscans .* O.nsess],'*.img',['Select images for Subject ' mySubject],pwd,0);
	end


	% ----------------------------------------------------------------------------------
	% Load/build the analysis matrix 
	% ----------------------------------------------------------------------------------
	disp(['	Building design matrices'])

	if desLoad == 1
		% ----------------------------------------------------------------------------------
		% Get Analysis file with design matrix
		% ----------------------------------------------------------------------------------

		if ~any(strcmp(N,'an_files')),O.an_files = [];,end

		if length(O.an_files) < O.snums(mySubIndex)
			O.an_files{O.snums(mySubIndex)} = [];
		end

		if isempty(O.an_files{O.snums(mySubIndex)})
			O.an_files{O.snums(mySubIndex)} = spm_get([1],'SPM.mat',['Select individual SPM.mat for Subject ' mySubject],pwd,0);
		end

		% ----------------------------------------------------------------------------------
		% Get SPM.mat file and build deconv matrix, if necessary
		% ----------------------------------------------------------------------------------

		eval(['load ' O.an_files{O.snums(mySubIndex)} ',''xX'',''Sess'';'])

		if dodeconv
			% concatenate stick functions for individual sessions
			c = [];		
			for kk = 1:length(Sess) 
				c = [c;Sess{kk}.sf];
			end
			c = cell2mat(c);
			O.evts_sf{O.snums(mySubIndex)} = c;

			if ~any(strcmp(N,'eres')), O.eres = input('Enter num of samples in sf per TR: ');,end
			if ~any(strcmp(N,'dxTRs')), O.dxTRs = input('Enter num of TRs to deconvolve: ');,end
			DX = tor_make_deconv_mtx2(O.evts_sf{O.snums(mySubIndex)},O.dxTRs,O.eres); 
		end

	elseif desLoad == 2
		% ----------------------------------------------------------------------------------
		% Get excel lists and turn them into a design
		% ----------------------------------------------------------------------------------
		warning('This is not implemented yet.  Should fit individual design matrices with SPM.')
		disp('...but deconvolution still done.')

		O = tor_build_sf(O.evt_onsets{mySubIndex},O);
		DX = O.DX;		% session-specific deconvolution matrix

	else disp('Unexpected desLoad value'), desLoad, return
	end



	% ----------------------------------------------------------------------------------
	% get the timeseries for each cluster
	% ----------------------------------------------------------------------------------
	fprintf(1,['	Getting timeseries data: Cluster '])

	sub_clusters = [];
	for myCl = 1:length(clusters)
		fprintf(1,[num2str(myCl) ' '])

		cl = clusters(myCl);


		O.coords = cl.XYZ';
		if ~isempty(O.im_files{O.snums(mySubIndex)}),
			try
				imn = O.im_files{O.snums(mySubIndex)};
				if ~iscell(imn)
					imn = mat2cell(imn,O.nscans * ones(O.nsess,1),size(imn,2)); 
 				end

				ts = timeseries2('multi',imn,O);
				cl.timeseries = ts.avg;
				cl.all_data = ts.indiv;
			catch
				lasterr
			end
		else
			disp('No timeseries data extracted - image names empty.')
		end

		sub_clusters = [sub_clusters, cl];

	end



	% ----------------------------------------------------------------------------------
	% Fit design matrix and deconvolve, if necessary
	% ----------------------------------------------------------------------------------
	fprintf(1,['	Analyzing timeseries data'])

	if desLoad == 1	% if you specified SPM.mat with xX
	   try
		% ----------------------------------------------------------------------------------
		% Analyze cluster ROIs with xX
		% ----------------------------------------------------------------------------------
		if ~any(strcmp(N,'TR')), O.TR = input('Enter TR: ');,end
		if ~isfield(xX,'RT'),xX.RT = O.TR;,end
		[sub_clusters] = analyze_cluster_rois(sub_clusters,xX,O.HP);
		for i = 1:length(sub_clusters), sub_clusters(i).xXb = sub_clusters(i).b;,end
	   catch
		disp('Error in analysis: skipping analysis.')
	   end

	end

	if dodeconv
		% for i = 1:length(sub_clusters), sub_clusters(i).xXb = sub_clusters(i).b;,end
		xDX.X = DX;
		xDX.K = ones(1,O.nsess);
	
		if ~any(strcmp(N,'TR')), O.TR = input('Enter TR: ');,end		

		xDX.RT = O.TR;
		[sub_clusters] = analyze_cluster_rois(sub_clusters,xDX,O.HP);
		for i = 1:length(sub_clusters), sub_clusters(i).DXb = sub_clusters(i).b;,end
	end
			

	GenericClusters = clusters;
	clusters = sub_clusters;
	eval(['save sub' mySubject '_clusters clusters xX DX'])
	clusters = GenericClusters;

end	% for each subject








catch
	disp('Problem in process batch_multisubject_cluster.')
	whos
	mySubIndex
	lasterr
	keyboard
end


return
	
	