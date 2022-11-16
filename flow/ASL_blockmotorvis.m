function ASL_blockmotorvis(root,msk,tr,t_off,t_on,t_start,doScl)

% Set default for root:
if nargin<1 || isempty(root)
    root = 'sub';
end

% Read root image from file
if contains(root,'.nii')
    [sub,fov,subtr] = nii2matrix(root);
elseif ~isempty(dir([root '*.nii*']))
    [sub,fov,subtr] = nii2matrix([root '.nii']);
else
    error('root: %s must point to a timeseries nifti file',root);
end

% Get dimensions
nframes = size(sub,4);
dim = [size(sub,1),size(sub,2),size(sub,3)];

% Set defaults
if nargin<2 || isempty(msk)
    % if msk is not defined, do not mask
    msk = ones(dim);
end
    
if nargin<3 || isempty(tr)
    tr = subtr;
end

if nargin<4 || isempty(t_off)
    t_off = 30;
end

if nargin<5 || isempty(t_on)
    t_on = t_off;
end

if nargin<6 || isempty(t_start)
    t_start = 0;
end

if nargin<7 || isempty(doScl)
    doScl = 1;
end

% Construct hi-res timing
t_total = tr*(nframes-1);
t_hires = 0:1e-3:t_total;
t_lowres = 0:tr:t_total;

% Calculate stimulus & hrf
stim = 1 * (t_hires>=t_start) .* (mod(t_hires-t_start,t_on+t_off) > t_off);
[~, p] = spm_hrf(tr*1e-3);
p(5) = 100;
[hrf, ~] = spm_hrf(tr*1e-3,p);

% Convolve stimulus with hrf
x_act_hires = conv(stim,hrf); x_act_hires = x_act_hires(1:length(t_hires));

% Interpolate hires stimulus to get regressors
x_act = interp1(t_hires,x_act_hires,t_lowres)';
x_flow = ones(nframes,1);

% Create design & contrast matrix
A = [x_act,x_flow];
Ncon = size(A,2);
C = eye(Ncon); % Estimate all contrasts individually

% Perform least-squares regression to find beta values
sub_cat = reshape(permute(sub,[4,1:3]),nframes,prod(dim));

beta_cat = pinv(A) * sub_cat; % beta estimates
t_cat = zeros(Ncon,prod(dim));
for c = 1:Ncon % loop through contrasts
    t_cat(c,:) = C(c,:) * beta_cat / std(C(c,:) * beta_cat);
end

% Extract contrast beta values
beta_act = reshape(beta_cat(1,:),dim).*msk;
t_act = reshape(t_cat(1,:),dim).*msk;
beta_flow = reshape(beta_cat(2,:),dim).*msk;
t_flow = reshape(t_cat(2,:),dim).*msk;

% Plot
figure
subplot(2,2,1:2)
plot(t_hires,x_act_hires,t_hires,ones(size(t_hires))); hold on
scatter(t_lowres,x_act); scatter(t_lowres,x_flow);
title('Regressors');
subplot(2,2,3)
lightbox(t_flow);
title('T-map Baseline flow');
subplot(2,2,4);
lightbox(t_act);
title('T-map Activation');

matrix2nii('./flow_betamap.nii',beta_flow,fov,tr,doScl);
matrix2nii('./flow_tmap.nii',t_flow,fov,tr,doScl);
matrix2nii('./activation_betamap.nii',beta_act,fov,tr,doScl);
matrix2nii('./activation_tmap.nii',t_act,fov,tr,doScl);

end