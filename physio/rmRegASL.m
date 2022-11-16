function rmRegASL (root, allDesMat, xclude)
% 
%   ---- DO NOT USE .  IT destroys the ASL signal!
% 
% rmRegASL (file_root, allslices_Desmats. ,[xclude])
%
% removes known trends in ASL time series of images.  Aimed at ASL Physio correction
% 
%   (c) 2007 Luis Hernandez-Garcia
%   University of Michigan
%   report bugs to:  hernan@umich.edu
%
% The program solves the linear model specified in DesMat for the
% parameters by ordinary least squares. It returns the residuals to a 4D image file
% called residuals.img
%
% This prorgram splits both the data and the regressors into odds and evens
% and fits them separately!
%
% xclude lets you estimate, but not remove regressors in the
% design matrix.  If ommited we default to 1 (not removing the first
% regressor, which is assumed to be the baseline)

if nargin==2
    xclude=1;
end

[all_data, h] = read_img_series(root);
Nframes = size(all_data,1);

% check to see if these are AVW or NIFTI files:
if ~isfield(h,'magic')
    h = avw2nii_hdr(h);
end

output = zeros(size(all_data));
fprintf('\n Done reading the data. crunching ...');
%hnames = dir(sprintf('%s*.hdr',root));
%h = read_hdr(hnames(1).name);
Spix = h.dim(2) *h.dim(3);



hv = h;
hv.dim(1) = 3;
hv.dim(5) = 1;
% compute the variance after removing regressors
var1 = (std(all_data,0,1)).^2;
write_nii('varBefore.nii',var1, hv,0);

% splitting data
cdata = all_data(1:2:end, :);
tdata = all_data(2:2:end, :);

fprintf('\nCleaning ODD time points ...') 
slice = 0;
pcount = inf;
for pix =1:size(cdata,2)

    if pcount > Spix
        pcount=1;
        slice = slice +1;
        % check to see whether we have a different design matrix for each
        % slice...
        if size(size(allDesMat),2) == 3
            DesMat = squeeze(allDesMat(slice,1:2:end,:));
        else
            DesMat = allDesMat(1:2:end,:);
        end
        % adjust the design matrix' size by taking only the end
        % i.e. - this clips off the physio data collected during disdaqs
        % DesMat = DesMat((end-Nframes+1):end,:);
        
        xtx_inv = pinv(DesMat);
        fprintf('\r  slice %d ...',slice);
    end

    pix_data = tdata(:,  pix );
    beta_est = xtx_inv*pix_data;

    %    if pix==5000
    %        keyboard
    %    end

    xDM = DesMat(:,xclude+1:end);
    xBeta = beta_est(xclude+1:end,:);
    output(1:2:end, pix ) = pix_data - xDM *xBeta ;
    pcount = pcount +1;

end

fprintf('\nCleaning EVEN time points ...')
slice = 0;
pcount = inf;
for pix =1:size(cdata,2)

    if pcount > Spix
        pcount=1;
        slice = slice +1;
        % check to see whether we have a different design matrix for each
        % slice...
        if size(size(allDesMat),2) == 3
            DesMat = squeeze(allDesMat(slice,2:2:end,:));
        else
            DesMat = allDesMat(2:2:end,:);
        end
        % adjust the design matrix' size by taking only the end
        % i.e. - this clips off the physio data collected during disdaqs
        % DesMat = DesMat((end-Nframes+1):end,:);
        
        xtx_inv = pinv(DesMat);
        fprintf('\r  slice %d ...',slice);
    end

    pix_data = tdata(:,  pix );
    beta_est = xtx_inv*pix_data;

    %    if pix==5000
    %        keyboard
    %    end

    xDM = DesMat(:,xclude+1:end);
    xBeta = beta_est(xclude+1:end,:);
    output(2:2:end, pix ) = pix_data - xDM *xBeta ;
    pcount = pcount +1;

end


% compute the variance after removing regressors
var2 = (std(output,0,1)).^2 ;
write_nii('varAfter.nii',var2,hv,0);

fprintf('\n Writing output residuals ...');
outh = h;
outh.dim(5) = size(all_data,1);


%write_hdr( 'residuals.hdr', outh);
%write_img( 'residuals.img', output, outh);
write_nii( 'residuals.nii', output, outh,0);

warning on

save regressors xDM DesMat allDesMat

fprintf('\n...Done');
return





