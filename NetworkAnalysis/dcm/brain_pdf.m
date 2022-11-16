function brain_pdf;


close all

xdim = 32;
ydim=32;
zdim=32;
tdim=200;
showMovie = 0;
ph = phantom3d(32);

Npix = xdim*ydim*zdim;

[xx yy zz] = ndgrid( ...
    linspace(-1, 1, xdim) , ...
    linspace(-1, 1, ydim) , ...
    linspace(-1, 1, zdim));

mygrid.xx = xx;
mygrid.yy = yy;
mygrid.zz = zz;


%% part one:  describe the two SPATIAL patterns of the PDF
% cluster 1:
xc = 0.5;
yc = -0.1;
zc = 0;
sigmax = 0.1;
sigmay = 0.1;
sigmaz = 0.1;

cluster1 = gaussianCluster( [xc yc zc] , [sigmax sigmay sigmaz], mygrid );

% cluster 2:
xc = -0.5;
yc = 0.5;
zc = 0;
sigmax = 0.1;
sigmay = 0.1;
sigmaz = 0.1;

cluster2 = gaussianCluster( [xc yc zc] , [sigmax sigmay sigmaz], mygrid );

% cluster 3:
xc = 0.1;
yc = 0.5;
zc = 0;
sigmax = 0.1;
sigmay = 0.05;
sigmaz = 0.05;

cluster3 = gaussianCluster( [xc yc zc] , [sigmax sigmay sigmaz], mygrid );

% cluster 4:
xc = -0.5;
yc = -0.5;
zc = 0;
sigmax = 0.1;
sigmay = 0.2;
sigmaz = 0.2;

cluster4 = gaussianCluster( [xc yc zc] , [sigmax sigmay sigmaz], mygrid );

% Ps is the spatial part of the Psi function
Ps{1} = cluster1 ;%+ cluster2;
Ps{2} = cluster3 ;%+ cluster4;

% Ps{1} = cluster2;
% Ps{2} = cluster4;

%{
lightbox(Ps{1});
figure
lightbox(Ps{2});
figure
%}
%% Part two: Now we generate the temporal patterns

% Here is a design matrix as a generative model
TR = 2;
exp_duration = tdim*TR;

onsets{1} = 1:30:exp_duration-30;
durations{1} = 10*ones(size(onsets{1}));

onsets{2} = 3:30:exp_duration - 30;
durations{2} = 5*ones(size(onsets{1}));

X = buildDesMat(TR, exp_duration, onsets, durations, 0);

rawdata = zeros( tdim,  xdim*ydim*zdim);

modelType=2;

switch modelType
    case 1
        
        Pt{1} = X(:,2);
        Pt{2} = X(:,3);
        
        %%
        for t=1:tdim
            rawdata(t,:) = ...
                Pt{1}(t) * reshape(Ps{1}, 1, xdim*ydim*zdim) + ...
                Pt{2}(t) * reshape(Ps{2}, 1, xdim*ydim*zdim);
            
            ov([], ...
                100*reshape(rawdata(t,:), xdim, ydim, zdim), ...
                round(xdim/2), round(ydim/2), round(zdim/2),...
                0); caxis([-0.1, 0.1]);
            
            %lightbox(  reshape(rawdata(t,:), xdim, ydim, zdim), [-1 1], 5 );
            drawnow
        end;
        
    case 2
        % the spatial patterns: Ps, serve as weighting factors
        % the desing matrix columns serve as external inputs
        
        % TR = 2;
        % dt = TR;
        % u = X(:,2);
        
        % Implement a slow Input process into node 1 only
        tdim = 500;
        u = zeros(tdim,1);
        u(51:200) = 1 - exp(-0.03*[0:149]);
        u(201:end) = u(200)*exp(-0.01*[0:299]);
        u = [u u]';
        dt = 1;
        
        
        Ps1 = reshape(Ps{1}, 1, xdim*ydim*zdim);
        Ps2 = reshape(Ps{2}, 1, xdim*ydim*zdim);
        
        Ps1 = Ps1/max(Ps1);
        Ps2 = Ps2/max(Ps2);
        
        % Coefficients for DCM
        % self influence terms
        A = diag([-0.095 -0.95]);
        
        % cross influence terms
        B = [0       0 ;
            0.5    -0.0  ] ;
        % input unfluence terms
        C = diag([-0.01  0]);
        
        rawdata = zeros( tdim,  xdim*ydim*zdim);
        % baseline values
        x = ones(2,tdim);
        x0  = [1 1]';  % baseline equilibrium value;
        rawdata(1,:) = 0.9; %Ps1+Ps2;
        
        for t=2:tdim

%            in this model the connection strenghts are weighted by the
%            spatial distribution of the model
%             dRaw_dt =  ...
%                 A(1) * Ps1 .*  rawdata(t-1,:)  + ... % self-influence loops
%                 A(2) * Ps2 .*  rawdata(t-1,:)  + ...
%                 B(1) * Ps2 * sum(rawdata(t-1,:) .* Ps1) + ...% cross terms :  modulation by other networks
%                 B(2) * Ps1 * sum(rawdata(t-1,:) .* Ps2) + ...
%                 C(1) * Ps1 *  u(t) ;% external input into network 1
%             
%             rawdata(t,:) = rawdata(t-1,:) + dRaw_dt*dt;
            
            % compute DCM for a simple model
            dx_dt =  ...
                A * (x(:,t-1) -x0) + ... % self-influence terms
                B * (x(:,t-1) -x0) + ... % cross influence terms
                C * u(:,t-1)  ; % external input into network 1
            
            x(:,t) = x(:, t-1) + dx_dt*dt;
            
            % apply DCM to the spatial clustering pattern.
            rawdata(t,:) = x(1,t)*Ps1 + x(2,t) * Ps2;
            
            % Display stuff:
            xyz = [round(xdim/2),  round(ydim/2),  round(zdim/2)];
            ind = sub2ind([xdim ydim zdim], xyz(1), xyz(2), xyz(3));
            
            if showMovie
                ov([], ...
                    ph + reshape(rawdata(t,:), xdim, ydim, zdim), ...
                    xyz(1), xyz(2), xyz(3),...
                    0); drawnow
            end
            
            
            drawnow
            
        end
        
        % each cluster's  time course:
        ind1 = find(Ps1>0.8);
        ind2 = find(Ps2>0.8);
        
        ts1 = mean(rawdata(:,ind1),2);
        ts2 = mean(rawdata(:,ind2),2);
      
        
        ov([], ...
            ph + reshape(rawdata(t,:), xdim, ydim, zdim), ...
            xyz(1), xyz(2), xyz(3),...
            0); 
        
        subplot(224)
        plot(ts1, 'b'); hold on;
        plot(ts2,'r');
        plot(u','k'); hold off;
        %legend('node 1', 'node 2', 'input into node 1',)
        
        % 
        figure
        ts1 = x(1,:);
        ts2 = x(2,:);
        t = linspace(0,21,length(ts1));
        subplot(211), plot(t, u); title('System input (e.g., drug)')
        subplot(212), plot(t,ts1); title('Node Activity')
        hold on
        subplot(212), plot(t, ts2,'r'); 
        legend ('Node 1 Baseline', 'Node 2 Baseline')
        hold off
        

end


return
%%


function X=buildDesMat(TR, exp_duration, onsets, durations, doASLmod)
% function X = buildDesMat(TR, exp_duration, onsets, durations, doASLmod)
%
%     (c) 2010 Luis Hernandez-Garcia @ UM
%     report bugs to hernan@umich.edu
%
%     this is a function to make design matrices, including the ASL
%     modulation if you want to analyze ASL data without subtraction
%
%     The baseline is the first column and gets added in for you.
%
%     TR :  scanner repetition time ...please give all data in seconds
%
%     exp_duration : total duration of the run
%
%     onsets: onset times for each block/event of the stimulation (or condition)
%         this should be a cell array containing the onset times of each condition like this:
%         onsets{1} = [ 1, 40, 60, 90 ....];
%         onsets{2} = [ 10, 20, ....];
%
%     durations : similarly, give the duration of each of the above blocks.  If they are events, make them 1
%         make sure that the sizes match up.
%         durations{1} = [20 20 20 ...];
%         durations{2} = [3 1 5 ....];
%
%     doASLmod : 1 if you're making a design matrix for an unssubtracted ASL experiment,
%                 0 otherwise
%

% Allocate space:
D = zeros(exp_duration , length(onsets)+1);
Nconds = length(onsets);

%% describe the basic set of conditions
D(:,1) = 1;
for cond=1:Nconds
    for blk = 1:length(onsets{cond})
        D( onsets{cond}(blk) : (onsets{cond}(blk) + durations{cond}(blk) ), cond+1) = 1;
    end
end

%% Now convolve with HRF:
h = spm_hrf(1);

for r=1:Nconds+1
    reg = D(:,r);
    reg = conv(reg,h);
    reg = reg(1:exp_duration);
    D(:,r) =  reg;
end
D(:,1) = 1;

%% resample to TR
tsamp1=1:exp_duration;
tsamp2=1:TR:exp_duration;
for r=1:Nconds+1
    reg = D(:,r);
    reg = interp1(tsamp1, reg, tsamp2);
    D2(:,r) =  reg;
end
D=D2;

%% now generate the ASL modulated version:
if doASLmod
    D2 = D;
    D2(1:2:end,:) = -D2(1:2:end,:);
    
    % concatenate the two:
    X = [D D2/2];
else
    X=D;
end

return



%%
function [fig1, fig2, fig3] =  ov(h,d,x,y,z,roi)
%function [fig1, fig2, fig3] =  ov(h,d,x,y,z,roi)
%
% (c) 2005 Luis Hernandez-Garcia
% University of Michigan
% report bugs to:  hernan@umich.edu
%

if isempty(h)
    fig1=subplot (223);
    imagesc(squeeze(d(:,:,z)))  ,axis xy %,axis tight
    
    hold on; plot(y,x,'go');hold off;
    set(gca,'Position',[ 0.05    0.01    0.45    0.45]);
    
    fig2=subplot (222);
    imagesc(squeeze(d(:,y,:))'),axis xy %, axis tight
    hold on; plot(x,z,'go');hold off;
    set(gca,'Position',[ 0.5    0.5    0.45    0.45]);
    
    fig3=subplot (221);
    imagesc(squeeze(d(x,:,:))'), axis xy%, axis tight
    hold on; plot(y,z,'go');hold off;
    set(gca,'Position',[ 0.05    0.5    0.45    0.45]);
    
    if roi>0
        value=mean(mean(d(x-roi:x+roi, y-roi:y+roi,z-roi:z+roi)));
    else
        value=d(x,y,z);
    end
    fprintf('\r(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, value);
else
    stretch = abs(h.zsize/h.xsize);
    %colordef black
    
    fig1=subplot (223);
    imagesc(squeeze(d(:,:,z))), axis ([0.5 h.ydim 0.5 h.xdim]) , axis xy %xy %,axis tight
    hold on; plot(y,x,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
    ht = abs(h.xdim*h.xsize);
    wt = abs(h.ydim*h.ysize);
    biggest = max([ht wt]);
    scl = 0.45 / biggest;
    set(gca,'Position',[ 0.01    0.05    wt*scl    ht*scl]);
    
    fig2=subplot (222);
    imagesc(squeeze(d(:,y,:))'), axis ([0.5 h.xdim 0.5 h.zdim]),axis  xy %, axis tight
    hold on; plot(x,z,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
    %set(gca,'Position',[ 0.5    0.5    0.45    0.45]);
    wt = abs(h.xdim*h.xsize);
    ht = abs(h.zdim*h.zsize);
    %biggest = max([ht wt]);
    scl = 0.45 / biggest;
    set(gca,'Position',[ 0.51    0.55    wt*scl    ht*scl]);
    
    fig3=subplot (221);
    imagesc(squeeze(d(x,:,:))'), axis ([0.5 h.ydim 0.5 h.zdim]), axis xy%, axis tight
    hold on; plot(y,z,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
    ht = abs(h.zdim*h.zsize);
    wt = abs(h.ydim*h.ysize);
    %biggest = max([ht wt]);
    scl = 0.45 / biggest;
    set(gca,'Position',[ 0.01    0.55    wt*scl    ht*scl]);
    
end

%

%str = sprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, tmp);
%subplot(221), title(str)

return


%%
function cluster = gaussianCluster( xyz, sigmas, mygrid)
% function cluster = gaussianCluster( xyz, sigmas, mygrid)
% This function constructs a Gaussian in 3D
%

xc = xyz(1);
yc = xyz(2);
zc = xyz(3);

sigmax = sigmas(1);
sigmay = sigmas(2);
sigmaz = sigmas(3);


cluster = 1/sqrt(2*pi*(sigmax^2 + sigmay^2 + sigmaz^2)) * exp(  ...
    -(mygrid.xx - xc ).^2 / (2*sigmax^2)  ...
    -(mygrid.yy - yc ).^2 / (2*sigmay^2)  ...
    -(mygrid.zz - zc ).^2 / (2*sigmaz^2) );

return
