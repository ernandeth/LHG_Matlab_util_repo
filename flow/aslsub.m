function aslsub(root,navg,first,last,order,sur,rel,tr,doScl)

% Set default for root:
if nargin<1 || isempty(root)
    root = 'timeseries_mag';
end

% Read root image from file
if contains(root,'.nii')
    [ts,fov,tstr] = nii2matrix(root);
elseif ~isempty(dir([root '*.nii*']))
    [ts,fov,tstr] = nii2matrix([root '.nii']);
else
    error('root: %s must point to a timeseries nifti file',root);
end

% Get dimensions
dim = [size(ts,1),size(ts,2),size(ts,3)];
nframes = size(ts,4);

% Set default for navg
if nargin<2 || isempty(navg)
    navg = 1;
end

% Set default for first
if nargin<3 || isempty(first)
    first = 3;
end

% Set default for last
if nargin<4 || isempty(last)
    last = nframes;
end

% Set default for order
if nargin<5 || isempty(order)
    order = 0;
end

% Set default for sur
if nargin<6 || isempty(sur)
    sur = 1;
end

% Set default for rel
if nargin<7 || isempty(rel)
    rel = 0;
end

% Set default for tr
if nargin<8 || isempty(tr)
    tr = tstr;
end

% Set default for doScl
if nargin<9 || isempty(doScl)
    doScl = 1;
end

% Split data into control and label
ts_con = ts(:,:,:,(first+1*~order):2:last);
ts_tag = ts(:,:,:,(first+1*order):2:last);

if ~sur
    % Subtract control/label (pairwise)
    ts_sub = ts_tag - ts_con;

    % Average together frames of sub
    if navg~=1 && mod(size(ts_sub,4),navg)==0
        sub_tmp = zeros([size(ts_sub,1),size(ts_sub,2),size(ts_sub,3),size(ts_sub,4)/navg]);
        ts_con_tmp = sub_tmp;
        ts_tag_tmp = sub_tmp;
        for n = 1:navg:size(ts_sub,4)
            sub_tmp(:,:,:,(n-1)/navg + 1) = mean(ts_sub(:,:,:,n:n-1+navg),4);
            ts_con_tmp(:,:,:,(n-1)/navg + 1) = mean(ts_con(:,:,:,n:n-1+navg),4);
            ts_tag_tmp(:,:,:,(n-1)/navg + 1) = mean(ts_tag(:,:,:,n:n-1+navg),4);
        end
        ts_sub = sub_tmp;
        ts_con = ts_con_tmp;
        ts_tag = ts_tag_tmp;
    elseif navg~=1
        warning('Cannot compute %d avgs for each point since number of timepoints in subtraction (%d) is not divisible by %d',...
            navg,nframes,navg);
    end
else
    % Subtract neighboring frame from each frame
    ts_sub = zeros([dim,last-first+1]);
    
    ts_sub(:,:,:,1) = (-1)^(order) * (ts(:,:,:,first) - ts(:,:,:,first+1));
    ts_sub(:,:,:,end) = (-1)^(last-first+order) * (ts(:,:,:,last) - ts(:,:,:,last-1));
    
    for tn = first+1:last-1
        imf = ts(:,:,:,tn);
        imnp = ts(:,:,:,tn-1);
        imnn = ts(:,:,:,tn+1);
        
        ts_sub(:,:,:,tn-first) = (-1)^(tn-first+order) * (2*imf - imnp - imnn);
    end
end

% If rel, divide subtraction timeseries by M0
if rel && first>1
    M0 = mean(ts(:,:,:,1:first-1),4) + eps();
    ts_sub = ts_sub./M0;
elseif rel
    warning('cannot compute subtractions relative to M0 without any M0 frames at beginning of timeseries');
end

% Calculate temporal mean and std of timeseries
mean_con = mean(ts_con,4); std_con = std(ts_con,[],4);
mean_tag = mean(ts_tag,4); std_tag = std(ts_tag,[],4);
mean_sub = mean(ts_sub,4); std_sub = std(ts_sub,[],4);

% Write images out to file
matrix2nii('./con.nii',ts_con,fov,tr*2,doScl);
matrix2nii('./tag.nii',ts_tag,fov,tr*2,doScl);
matrix2nii('./sub.nii',ts_sub,fov,sur*tr + 2*~sur*tr,doScl); % sur mode conserves tr

matrix2nii('./mean_con.nii',mean_con,fov,tr*2,doScl);
matrix2nii('./mean_tag.nii',mean_tag,fov,tr*2,doScl);
matrix2nii('./mean_sub.nii',mean_sub,fov,sur*tr + 2*~sur*tr,doScl); % sur mode conserves tr

matrix2nii('./std_con.nii',std_con,fov,tr*2,doScl);
matrix2nii('./std_tag.nii',std_tag,fov,tr*2,doScl);
matrix2nii('./std_sub.nii',std_sub,fov,sur*tr + 2*~sur*tr,doScl); % sur mode conserves tr

end
