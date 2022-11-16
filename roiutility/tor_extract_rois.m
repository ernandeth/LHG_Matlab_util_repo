function [clusters,SPM,xX,xCon] = tor_extract_rois(imnames,varargin)
% function [clusters, SPM, xX, xCon] = tor_extract_rois(imnames [can be empty],[opt] SPM, [opt] VOL, [opt] xX)
%
% this function gets timeseries data from all clusters in an SPM results output.
% input: 
%	imnames: a matrix of image names, in spm_list_files output format
%		 if empty, no timeseries data will be extracted.
%	SPM:	 SPM variable from loaded results
%	VOL:	 VOL variable from loaded results
%	[Last 2 arguments are optional.  Use if results are already loaded into workspace]
%
% Automatic fitting of model to cluster timeseries average using analyze_cluster_rois
% with High-Pass filter length of your choice.
% This only works if you input only the file names or input all optional arguments, including xX
%
% 10/17/01 by Tor Wager


	if nargin == 1
		% ----------------------------------------------------------------------------------
		% get SPM info from SPM.mat
		% ----------------------------------------------------------------------------------
		[SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
	elseif nargin == 3
		SPM = varargin{1};
		VOL = varargin{2};
	elseif nargin == 4
		SPM = varargin{1};
		VOL = varargin{2};
		xX = varargin{3};
	elseif nargin ~= 3
		error('Wrong number of arguments.  Use 1 if loading from SPM.mat, 3 args for no analysis, 4 with analysis.')
	end
		
	% ----------------------------------------------------------------------------------
	% get cluster index number from each voxel
	% ----------------------------------------------------------------------------------
	cl_index = spm_clusters(SPM.XYZ);
	
	% ----------------------------------------------------------------------------------
	% define each cluster as cell in array.
	% ----------------------------------------------------------------------------------
	clusters = [];
	for i = 1:max(cl_index)
		disp(['Extracting cluster ' num2str(i) ' of ' num2str(max(cl_index))])
		a = find(cl_index == i);

		cl.title = SPM.title;
		cl.threshold = SPM.u;
		cl.voxSize = VOL.VOX;
		cl.name = [cl.title '_' num2str(i) '_' mat2str(size(a,2)) '_voxels'];
		cl.numVox = size(a,2);
		cl.Z = SPM.Z(a);
		cl.XYZmm = SPM.XYZmm(:,a);
		cl.XYZ = SPM.XYZ(:,a);
		try
			cl.pVoxelLev = spm_P(1,0,max(cl.Z),SPM.df,SPM.STAT,VOL.R,SPM.n);
			cl.pClustLev = spm_P(1,cl.numVox/prod(VOL.FWHM),SPM.u,SPM.df,SPM.STAT,VOL.R,SPM.n);
		catch
			warning('Can''t get SPM p voxel and cluster values.  Skipping spm_P.')
		end

		% ----------------------------------------------------------------------------------
		% get the timeseries for the cluster
		% ----------------------------------------------------------------------------------
		O.coords = cl.XYZ';
		if ~isempty(imnames),
			try
				ts = timeseries2('multi',imnames,O);
				cl.timeseries = ts.avg;
				cl.all_data = ts.indiv;
			catch
			end
		else
			disp('No timeseries data extracted - image names empty.')
		end


		clusters = [clusters, cl];
	end

	if ~isempty(imnames) & exist('xX') == 1

		% ----------------------------------------------------------------------------------
		% adjust timeseries for each cluster and fit xX.X model to timeseries data
		% ----------------------------------------------------------------------------------
		try
			clusters = analyze_cluster_rois(clusters,xX);
		catch
			disp(['Error analyzing timeseries clusters - skipping analysis.'])
		end

	end
		
	% ----------------------------------------------------------------------------------
	% save this data in current directory
	% ----------------------------------------------------------------------------------
	matname = ['clusters_' deblank(SPM.title(~(SPM.title==' ')))];
	str = ['save ' matname ' clusters'];
	eval(str)



return

